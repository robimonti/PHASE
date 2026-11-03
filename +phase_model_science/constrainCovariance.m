function [coefficients,info]=constrainCovariance(coefficients,fun,lag,empirical,cofactor,model,options,direction)
%CONSTRAINCOVARIANCE Optional bounded WLS refinement of the native kernel.
% Range = e-folding length (Gaussian/cosine envelope), not first zero crossing.
options=phase_model_science.options(options);
prior=options.(['correlation' direction]); maximum=options.(['correlation' direction 'Max']);
info=struct('enabled',false,'definition','e-folding length of the covariance envelope');
if isnan(prior)&&isnan(maximum),return;end
power=2;if strcmp(model,'exponential'),power=1;end
initial=coefficients(:);nativeRange=initial(2)^(-1/power);
lo=eps;hi=Inf;
if isfinite(prior)
    fraction=options.correlationTolerance/100;
    lo=prior*(1-fraction);hi=prior*(1+fraction);
end
if isfinite(maximum),hi=min(hi,maximum);end
if isnan(prior),prior=min(nativeRange,hi);end
initial(2)=max(lo,min(hi,prior))^(-power);
% Lag zero contains the nugget, which is estimated separately by STmodel.
valid=isfinite(lag)&isfinite(empirical)&isfinite(cofactor)&cofactor>0&lag>0;
if nnz(valid)<3,error('PHASE:correlationRange','At least three nonzero covariance bins are needed for manual range fitting.');end
x=lag(valid);y=empirical(valid);weight=1./sqrt(cofactor(valid));weight=weight/max(weight);
scale=max(abs(initial),eps);start=initial./scale;
lower=zeros(size(initial));upper=inf(size(initial));
lower(1)=eps;lower(2)=hi^(-power);upper(2)=lo^(-power);
% Keep the native cosine frequency: this control changes range, not period.
if numel(initial)>2,lower(3:end)=initial(3:end);upper(3:end)=initial(3:end);end
objective=@(p) weight.*(fun(p.*scale,x)-y);
[fitted,~,~,flag]=lsqnonlin(objective,start,lower./scale,upper./scale, ...
    optimoptions('lsqnonlin','Display','off','MaxIterations',500));
if flag<=0,error('PHASE:correlationRange','Bounded covariance fit did not converge. Review the supplied range.');end
coefficients=fitted.*scale;
info.enabled=true;info.nativeRange=nativeRange;info.suggestedRange=prior;
info.minimumRange=lo;info.maximumRange=hi;info.fittedRange=coefficients(2)^(-1/power);
info.atBound=abs(info.fittedRange-lo)<=1e-4*max(1,lo)||abs(info.fittedRange-hi)<=1e-4*max(1,hi);
end
