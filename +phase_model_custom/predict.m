function prediction = predict(q,referenceDay,c,ids,fit,residual)
%PREDICT Two-stage estimate and covariance, including shared-data terms.
% F=Jq*L, G=W*(I-J*L), total influence H=F+G.
% STC error covariance: Kqq - H*Koq - Kqo*H' + H*Cyy*H'.
% DET sampling covariance: H*N*H'. Nonlinear results use local derivatives.
reference=q; reference(:,4)=referenceDay;
[m,Jall]=phase_model_custom.evaluate(fit.definition,fit.values,q,fit.context);
[m0,J0all]=phase_model_custom.evaluate(fit.definition,fit.values,reference,fit.context);
Jq=Jall(:,fit.free)-J0all(:,fit.free);
m=m-m0;
if strcmp(residual.kind,'STC')
    K=phase_model_custom.kernel(q,c,residual.settings) ...
        -phase_model_custom.kernel(reference,c,residual.settings);
    W=(residual.factor'\(residual.factor\K'))';
    stochastic=K*residual.alpha;
    t=q(:,4); t0=referenceDay;
    if isfield(residual,'native') && residual.native
        zero=residual.settings.timeFunction(residual.settings.timeParameters,0)* ...
            residual.settings.spaceFunction(residual.settings.spaceParameters,0)/residual.settings.timeParameters(1);
        temporal=residual.settings.timeFunction(residual.settings.timeParameters,abs(t-t0));
        targetVariance=2*(zero-temporal*residual.settings.spaceFunction(residual.settings.spaceParameters,0)/residual.settings.timeParameters(1));
    else
        d=abs(t-t0)/residual.settings.timeRange;
        if strcmp(residual.settings.timeKernel,'gaussian'), corr=exp(-d.^2); else, corr=exp(-d); end
        targetVariance=2*residual.settings.signalVariance*(1-corr);
    end
else
    B=residual.basis(q,residual.nativeState)-residual.basis(reference,residual.nativeState);
    W=B*residual.coefficientMap;
    stochastic=B*residual.coefficients;
    K=zeros(size(W)); targetVariance=zeros(size(q,1),1);
end
F=Jq*fit.L;
G=W-(W*fit.J)*fit.L;
H=F+G;
if strcmp(residual.kind,'STC')
    HC=H*residual.C; FC=F*residual.C; GC=G*residual.C;
else
    HC=H*residual.noiseVariance;
    FC=F*residual.noiseVariance;
    GC=G*residual.noiseVariance;
end
total=targetVariance+sum(HC.*H,2)-2*sum(H.*K,2);
deterministicVariance=sum(FC.*F,2);
residualVariance=targetVariance+sum(GC.*G,2)-2*sum(G.*K,2);
cross=sum(FC.*G,2)-sum(F.*K,2);
if isfield(fit.definition,'assignedCovariance')
    P=fit.definition.assignedCovariance;
    [~,Jobs]=phase_model_custom.evaluate(fit.definition,fit.values,c,fit.context);
    sensitivity=zeros(numel(fit.free),numel(fit.values));
    if isfield(fit,'assignedSensitivity'),sensitivity=fit.assignedSensitivity;end
    D=Jall-J0all+Jq*sensitivity;
    R=-W*(Jobs+fit.J*sensitivity); Q=D+R;
    total=total+sum((Q*P).*Q,2);
    deterministicVariance=deterministicVariance+sum((D*P).*D,2);
    residualVariance=residualVariance+sum((R*P).*R,2);
    cross=cross+sum((D*P).*R,2);
end
scale=max([1;abs(targetVariance);abs(sum(HC.*H,2))]);
if any(total < -1e-7*scale)
    error('PHASE:customVariance','Negative propagated variance beyond numerical precision.');
end
prediction=struct('total',m+stochastic,'deterministic',m,'residual',stochastic, ...
    'variance',max(total,0),'deterministicVariance',max(deterministicVariance,0), ...
    'residualVariance',residualVariance,'crossCovariance',cross);
end
