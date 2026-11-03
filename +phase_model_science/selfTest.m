function selfTest()
%SELFTEST Independent algebraic checks for trend removal/restoration and error.
addpath(fullfile(phase_model_beta.projectRoot(),'MatlabFunctions'));
c=phase_model_beta.defaultConfig();c.temporalTrend='quadratic';
c.customExpression='a*t+b.*sin(t/90)';
c.customParameters='[{"name":"a","value":0.02,"estimate":true,"std":null},{"name":"b","value":0.2,"estimate":false,"std":0.1}]';
d=phase_model_custom.parseDefinition(c);times=0:30:600;
[s,t]=ndgrid([0 40 90],times);points=[s(:),zeros(numel(s),1),s(:),t(:)];
context=struct('origin',[0 0],'referenceDay',180,'times',times, ...
    'xy',[0 0;40 0;90 0]);ids=repmat((1:3)',numel(times),1);
y=0.03*(points(:,4)-180)+0.0002*(points(:,4).^2-180^2) ...
    +0.2*(sin(points(:,4)/90)-sin(2));
fit=phase_model_custom.fitDeterministic(d,points,y,context);
fit=phase_model_custom.augmentTemporalTrend(fit,points,y,ids,c,'');
assert(max(abs(y-fit.mean))<1e-6,'Known deterministic trend was not removed.');
assert(fit.temporalTrend.degree==2);
% Numerical derivative of BOTH stages, including their shared observations.
changed=y;changed(5)=changed(5)+1e-4;
refit=phase_model_custom.fitDeterministic(d,points,changed,context);
refit=phase_model_custom.augmentTemporalTrend(refit,points,changed,ids,c,'');
actual=[refit.values(refit.free)-fit.values(fit.free); ...
    refit.temporalTrend.coefficients-fit.temporalTrend.coefficients]/1e-4;
assert(max(abs(actual-fit.L(:,5)))<1e-3);
settings=struct('dimension','1D','native',true,'timeParameters',[2,1/70^2], ...
    'spaceParameters',[2,1/80^2],'timeFunction',@(p,h)p(1)*exp(-p(2)*h.^2), ...
    'spaceFunction',@(p,h)p(1)*exp(-p(2)*h.^2));
C=phase_model_custom.kernel(points,points,settings)+0.2*eye(numel(y));
residual=struct('kind','STC','native',true,'settings',settings,'C',C, ...
    'factor',chol(C,'lower'),'alpha',C\(y-fit.mean));
q=[20 0 20 550;70 0 70 90;20 0 20 0];q0=q;q0(:,4)=0;
p=phase_model_custom.predict(q,0,points,ids,fit,residual);
expected=0.03*q(:,4)+0.0002*q(:,4).^2+0.2*sin(q(:,4)/90);
assert(max(abs(p.total-expected))<1e-6,'Trend was not restored at the requested reference.');
assert(max(abs(p.total-p.customDeterministic-p.temporalTrend-p.residual))<1e-10);
[~,J]=phase_model_custom.evaluate(d,fit.values,q,context);
[~,J0]=phase_model_custom.evaluate(d,fit.values,q0,context);
Jq=[J(:,fit.free)-J0(:,fit.free),phase_model_science.trendBasis(q(:,4),0,fit.temporalTrend.scale,2)];
K=phase_model_custom.kernel(q,points,settings)-phase_model_custom.kernel(q0,points,settings);
W=K/C;H=W+(Jq-W*fit.J)*fit.L;
Kqq=phase_model_custom.kernel(q,q,settings)-phase_model_custom.kernel(q,q0,settings) ...
    -phase_model_custom.kernel(q0,q,settings)+phase_model_custom.kernel(q0,q0,settings);
[~,Jobs]=phase_model_custom.evaluate(d,fit.values,points,context);
Q=J-J0+Jq*fit.assignedSensitivity-W*(Jobs+fit.J*fit.assignedSensitivity);
oracle=diag(Kqq-H*K'-K*H'+H*C*H'+Q*d.assignedCovariance*Q');
assert(max(abs(p.variance-oracle))<1e-8,'Trend propagation lost a cross-covariance term.');
assert(max(abs(p.variance-p.deterministicVariance-p.residualVariance-2*p.crossCovariance))<1e-8);
assert(p.total(3)==0&&p.variance(3)==0);
% Blocked validation recognizes a large quadratic mean; a zero series is unchanged.
info=phase_model_science.temporalDiagnostic(times,repmat(0.0002*times.^2,3,1));
assert(info.suggestedDegree==2);
info=phase_model_science.temporalDiagnostic(times,zeros(3,numel(times)));
assert(info.suggestedDegree==0);
% Bounds use e-folding lengths for all native kernel families.
lags=(0:10:300)';weights=ones(size(lags));
for family={'gaussian','exponential','gaussian_cosine'}
    model=family{1};power=2;
    if strcmp(model,'exponential'),power=1;fun=@(p,h)p(1)*exp(-p(2)*h);
    elseif strcmp(model,'gaussian'),fun=@(p,h)p(1)*exp(-p(2)*h.^2);
    else,fun=@(p,h)p(1)*exp(-p(2)*h.^2).*cos(p(3)*h*pi/2);end
    coeff=[3;1/150^power];if strcmp(model,'gaussian_cosine'),coeff(3)=0.01;end
    [native,info]=phase_model_science.constrainCovariance(coeff,fun,lags,fun(coeff,lags),weights,model,struct(),'Space');
    assert(isequal(native,coeff)&&~info.enabled);
    opts=struct('correlationSpace',50,'correlationTolerance',20,'correlationSpaceMax',55);
    [bounded,info]=phase_model_science.constrainCovariance(coeff,fun,lags,fun(coeff,lags),weights,model,opts,'Space');
    assert(info.enabled&&bounded(2)^(-1/power)>=40-1e-6&&bounded(2)^(-1/power)<=55+1e-6);
end
disp('Residual controls passed: remove/restore, shared-data uncertainty, blocked diagnostics, native/manual ranges.');
end
