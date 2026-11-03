function residualFieldSelfTest()
%RESIDUALFIELDSELFTEST Spatial generalisation, restoration and covariance oracle.
root=phase_model_beta.projectRoot();addpath(fullfile(root,'MatlabFunctions'));
[x,y]=ndgrid(linspace(-100,100,6),linspace(-70,70,6));
xy=[x(:) y(:)];times=0:30:540;np=size(xy,1);nt=numel(times);
c=[repmat(xy,nt,1),repmat((1:np)',nt,1),repelem(times(:),np)];
ids=repmat((1:np)',nt,1);ref=240;
% A field invisible to the AOI mean, but predictable at unseen positions.
r=.0002*c(:,1).*(c(:,4)-ref);
[state,B,M]=phase_model_science.selectResidualField(c,r,ids,'2D',ref);
assert(state.degree==1&&state.spatialDegree==1,'Did not choose the simplest exact spatial trend.');
assert(max(abs(B*state.coefficients-r))<1e-8);
assert(max(abs(M*r-state.coefficients))<1e-10);
assert(state.selected.timeRatio<.1&&state.selected.spaceRatio<.1);
q=[30 15 0 500;30 15 0 0];
Bq=phase_model_science.spaceTimeBasis(q,state,0);
assert(max(abs(Bq*state.coefficients-.0002*q(:,1).*q(:,4)))<1e-8);
% A native local field must also be predictable/restorable, not just a plane.
local=state;local.spatialKind='native cubic splines';local.spatialIntervals=2;
Alocal=phase_model_science.spaceTimeBasis(c,local,ref);
localCoefficients=50*sin((1:size(Alocal,2))');
localSignal=Alocal*localCoefficients;
[native,Anative,~]=phase_model_science.selectResidualField(c,localSignal,ids,'2D',ref);
assert(native.degree>0&&strcmp(native.spatialKind,'native cubic splines'),native.reason);
assert(max(abs(Anative*native.coefficients-localSignal))<1e-6);
[none,~,~]=phase_model_science.selectResidualField(c,zeros(size(r)),ids,'2D',ref);
assert(none.degree==0);
% Full two-stage influence, including an estimated upstream custom parameter.
config=phase_model_beta.defaultConfig();config.projDim='2D';
config.customExpression='a*t+b.*sin(t/100)';
config.customParameters='[{"name":"a","value":0,"estimate":true,"std":null},{"name":"b","value":0.2,"estimate":false,"std":0.1}]';
definition=phase_model_custom.parseDefinition(config);
context=struct('origin',[0 0],'referenceDay',ref,'times',times,'xy',xy);
yobs=.02*(c(:,4)-ref)+r+.2*(sin(c(:,4)/100)-sin(ref/100));
fit=phase_model_custom.fitDeterministic(definition,c,yobs,context);
fit=phase_model_custom.augmentTemporalTrend(fit,c,yobs,ids,config,'');
assert(fit.temporalTrend.spatialDegree==1);
assert(max(abs(yobs-fit.mean))<1e-7);
perturbed=yobs;perturbed(40)=perturbed(40)+1e-4;
changed=phase_model_custom.fitDeterministic(definition,c,perturbed,context);
changed=phase_model_custom.augmentTemporalTrend(changed,c,perturbed,ids,config,'');
actual=[changed.values(changed.free)-fit.values(fit.free);changed.temporalTrend.coefficients-fit.temporalTrend.coefficients]/1e-4;
assert(max(abs(actual-fit.L(:,40)))<1e-4,'The derivative of the two-stage field fit is incorrect.');
settings=struct('dimension','2D','native',true,'timeParameters',[2,1/70^2], ...
    'spaceParameters',[2,1/80^2],'timeFunction',@(p,h)p(1)*exp(-p(2)*h.^2), ...
    'spaceFunction',@(p,h)p(1)*exp(-p(2)*h.^2));
C=phase_model_custom.kernel(c,c,settings)+.2*eye(size(c,1));
residual=struct('kind','STC','native',true,'settings',settings,'C',C, ...
    'factor',chol(C,'lower'),'alpha',C\(yobs-fit.mean));
% Both targets must obey the full covariance oracle, including assigned errors.
for restore=[true false]
fit.temporalTrend.restore=restore;
p=phase_model_custom.predict(q,0,c,ids,fit,residual);
expected=fit.values(fit.free)*500+double(restore)*.0002*30*500+.2*sin(5);
assert(abs(p.total(1)-expected)<1e-7&&p.total(2)==0);
assert(max(abs(p.total-p.customDeterministic-p.restoredTrend-p.residual))<1e-10);
q0=q;q0(:,4)=0;
Jq=[q(:,4),double(restore)*phase_model_science.spaceTimeBasis(q,fit.temporalTrend,0)];
K=phase_model_custom.kernel(q,c,settings)-phase_model_custom.kernel(q0,c,settings);
H=K/C+(Jq-(K/C)*fit.J)*fit.L;
Kqq=phase_model_custom.kernel(q,q,settings)-phase_model_custom.kernel(q,q0,settings) ...
    -phase_model_custom.kernel(q0,q,settings)+phase_model_custom.kernel(q0,q0,settings);
[~,Jall]=phase_model_custom.evaluate(definition,fit.values,q,context);
[~,J0all]=phase_model_custom.evaluate(definition,fit.values,q0,context);
[~,Jobs]=phase_model_custom.evaluate(definition,fit.values,c,context);
Q=Jall-J0all+Jq*fit.assignedSensitivity-(K/C)*(Jobs+fit.J*fit.assignedSensitivity);
oracle=diag(Kqq-H*K'-K*H'+H*C*H'+Q*definition.assignedCovariance*Q');
assert(max(abs(p.variance-oracle))<1e-8);
assert(max(abs(p.variance-p.deterministicVariance-p.residualVariance-2*p.crossCovariance))<1e-8);
% DET uses the same target but a sampling covariance, tested independently.
map=pinv([ones(size(c,1),1),c(:,4)]);
det=struct('kind','DET','basis',@(v,state)[ones(size(v,1),1),v(:,4)], ...
    'nativeState',struct(),'coefficientMap',map,'coefficients',map*(yobs-fit.mean),'noiseVariance',.2);
pdet=phase_model_custom.predict(q,0,c,ids,fit,det);
W=([ones(size(q,1),1),q(:,4)]-[ones(size(q,1),1),q0(:,4)])*map;
Hd=Jq*fit.L+W-(W*fit.J)*fit.L;
Qd=Jall-J0all+Jq*fit.assignedSensitivity-W*(Jobs+fit.J*fit.assignedSensitivity);
oracleD=diag(.2*(Hd*Hd')+Qd*definition.assignedCovariance*Qd');
assert(max(abs(pdet.variance-oracleD))<1e-8);
end
% Automatic native check must detect local drifts hidden in the AOI mean.
localResidual=reshape(r,np,nt);
info=phase_model_science.temporalDiagnostic(times,localResidual,'');
assert(info.suggestedDegree==0&&info.suggestedPointDegree==1);
[retry,protected,degree]=phase_model_science.nativeTrendCheck(times,localResidual, ...
    phase_model_science.options(struct()),0,0,'','');
assert(retry&&protected==1&&degree==2);
% Degenerate geometries must not admit dependent added interactions.
u=repelem((0:3)',2);v=repmat([-1;1],4,1);
A=[ones(8,1),v,u,v.^2,u.^2,v.*u];
mask=phase_model_science.polynomialMask(6,0,1,'1D',A);
assert(rank(A(:,mask))==nnz(mask)&&mask(3)&&mask(6));
% 1D uses chainage, not east/north coordinates.
line=c;line(:,3)=c(:,1);
[one,A,~]=phase_model_science.selectResidualField(line,r,ids,'1D',ref);
assert(one.degree==1&&one.spatialDegree==1&&max(abs(A*one.coefficients-r))<1e-7);
% Do not manufacture a deterministic correction for independent noise.
stream=RandStream('mt19937ar','Seed',715);
noise=20*randn(stream,size(c,1),1);
[quiet,~,~]=phase_model_science.selectResidualField(c,noise,ids,'2D',ref);
assert(quiet.degree==0,'Independent noise was incorrectly assigned a deterministic field.');
% A resolved multi-cycle signal must not be flattened by the polynomial mean.
oscillatory=8*(1+.003*c(:,1)).*sin(2*pi*(c(:,4)-ref)/120);
[cycle,Acycle,~]=phase_model_science.selectResidualField(c,oscillatory,ids,'2D',ref);
if cycle.degree>0
    retained=oscillatory-Acycle*cycle.coefficients;
    assert(norm(retained)>.9*norm(oscillatory),'An oscillatory residual was absorbed into the mean.');
end
% A valid covariance can decay monotonically throughout a finite observation
% window: long correlation length does not imply a deterministic trend.
longRange=exp(-abs(times'-times)/(10*range(times)));
assert(all(diff(longRange(1,:))<0)&&longRange(1,end)>.9);
[~,bad]=chol(longRange);assert(bad==0);
disp('Automatic residual field passed: 1D/2D selection, held-out prediction, reference, derivatives and covariance oracle.');
end
