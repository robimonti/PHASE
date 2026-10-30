function result=scientificSelfTest()
%SCIENTIFICSELFTEST Local propagation oracle, assigned uncertainty and mean.
root=phase_model_beta.projectRoot();addpath(fullfile(root,'MatlabFunctions'));
stream=RandStream('mt19937ar','Seed',7101);
config=phase_model_beta.defaultConfig();config.customParameters='[{"name":"a","value":0.02,"estimate":true}]';
definition=phase_model_custom.parseDefinition(config);
context=struct('origin',[0,0],'referenceDay',30);
[s,t]=ndgrid([0,40,90],[0,15,45,60]); c=[s(:),zeros(numel(s),1),s(:),t(:)];
ids=repmat((1:3)',4,1); y=0.03*(c(:,4)-30);
fit=phase_model_custom.fitDeterministic(definition,c,y,context);
assert(abs(fit.values-0.03)<1e-6);
settings=struct('dimension','1D','native',true,'timeParameters',[2,1/35^2], ...
 'spaceParameters',[2,1/80^2],'timeFunction',@(p,d)p(1)*exp(-p(2)*d.^2), ...
 'spaceFunction',@(p,d)p(1)*exp(-p(2)*d.^2));
C=phase_model_custom.kernel(c,c,settings)+0.2*eye(numel(y)); L=chol(C,'lower');
residual=struct('kind','STC','native',true,'settings',settings,'C',C,'factor',L,'alpha',C\(y-fit.mean));
q=[20,0,20,55;70,0,70,10;20,0,20,0]; q0=q;q0(:,4)=0;
p=phase_model_custom.predict(q,0,c,ids,fit,residual);
K=phase_model_custom.kernel(q,c,settings)-phase_model_custom.kernel(q0,c,settings);
W=K/C; H=q(:,4)*fit.L+W*(eye(numel(y))-fit.J*fit.L);
Kqq=phase_model_custom.kernel(q,q,settings)-phase_model_custom.kernel(q,q0,settings) ...
-phase_model_custom.kernel(q0,q,settings)+phase_model_custom.kernel(q0,q0,settings);
oracle=diag(Kqq-H*K'-K*H'+H*C*H');
assert(max(abs(p.variance-oracle))<1e-8);
assert(max(abs(p.variance-p.deterministicVariance-p.residualVariance-2*p.crossCovariance))<1e-8);
joint=[C,K';K,Kqq];[V,D]=eig((joint+joint')/2);
draws=V*diag(sqrt(max(diag(D),0)))*randn(stream,size(joint,1),20000);
errors=H*draws(1:numel(y),:)-draws(numel(y)+1:end,:);
assert(max(abs(var(errors,0,2)-oracle)./max(oracle,0.05))<0.04);
assert(p.variance(3)==0&&p.total(3)==0);
% External assigned parameter errors include their effect on residual fitting.
config.customParameters='[{"name":"a","value":0.03,"estimate":false,"std":0.02}]';
fixed=phase_model_custom.fitDeterministic(phase_model_custom.parseDefinition(config),c,y,context);
pf=phase_model_custom.predict(q,0,c,ids,fixed,residual);
expected=diag(Kqq-W*K'-K*W'+W*C*W')+0.02^2*(q(:,4)-W*(c(:,4)-30)).^2;
assert(max(abs(pf.variance-expected))<1e-8);
assert(max(abs(pf.variance-pf.deterministicVariance-pf.residualVariance-2*pf.crossCovariance))<1e-8);
config.customResidualMean='constant';
mixed=phase_model_custom.fitDeterministic(phase_model_custom.parseDefinition(config),c,y,context);
pm=phase_model_custom.predict(q,0,c,ids,mixed,residual);
Jexternal=c(:,4)-30;
Hm=W*(eye(numel(y))-ones(numel(y))/numel(y));
expected=diag(Kqq-Hm*K'-K*Hm'+Hm*C*Hm')+ ...
    0.02^2*(q(:,4)-W*(Jexternal-mean(Jexternal))).^2;
assert(max(abs(pm.variance-expected))<1e-7,'Mixed assigned/estimated parameter propagation failed.');
assert(abs(mixed.assignedSensitivity(1,1)+mean(Jexternal))<1e-6);
% Optional constant is fitted jointly, never secretly removed by covariance.
config.customParameters='[{"name":"a","value":0.02,"estimate":true}]';
config.customResidualMean='constant';
meanFit=phase_model_custom.fitDeterministic(phase_model_custom.parseDefinition(config),c,y+7,context);
assert(abs(mean(meanFit.mean-y)-7)<1e-5&&abs(mean(y+7-meanFit.mean))<1e-6);
assert(abs(meanFit.values(1)-0.03)<1e-6);
% Reject incompatible assigned covariance and estimated parameters.
config.customParameterCovariance='[[1]]'; mustFail(@()phase_model_custom.parseDefinition(config));
config.customParameterCovariance='[]';config.customResidualMean='model';
config.customExpression='a.*(-expm1(-t./tau))';
config.customParameters='[{"name":"a","value":3,"estimate":true,"lower":0,"upper":20},{"name":"tau","value":20,"estimate":true,"lower":1,"upper":150}]';
d=phase_model_custom.parseDefinition(config);yn=phase_model_custom.evaluate(d,[5;35],c,context);
nonlinear=phase_model_custom.fitDeterministic(d,c,yn,context);
assert(max(abs(nonlinear.values-[5;35]))<1e-3);
changed=yn;changed(2)=changed(2)+1e-4;
refit=phase_model_custom.fitDeterministic(d,c,changed,context);
assert(max(abs((refit.values-nonlinear.values)/1e-4-nonlinear.L(:,2)))<0.02);
% Supplied DSC example, using its printed coefficients and original units.
addpath(fullfile(root,'examples','custom_models'));
assert(abs(DSC_quadratic_example(-0.01884026,-0.00542769,0,540,struct())-127.96285814570561)<1e-8);
config.customFunctionFile=fullfile(root,'examples','custom_models','DSC_quadratic_example.m');
config.customParameters='[]';d=phase_model_custom.parseDefinition(config);
v=phase_model_custom.evaluate(d,[],[-0.01884026,-0.00542769,0,540],struct('origin',[0,0],'referenceDay',0));
assert(abs(v-127.96313029450918)<1e-8);
result=struct('ok',true,'message','Scientific tests passed: covariance oracle, Monte Carlo, assigned uncertainties, residual mean, nonlinear fit and DSC example.');
disp(result.message);
end
function mustFail(action)
failed=false;try,action();catch,failed=true;end
assert(failed,'Invalid scientific model accepted.');
end
