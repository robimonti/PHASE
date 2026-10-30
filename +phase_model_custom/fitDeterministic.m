function fit = fitDeterministic(definition,coordinates,observations,context)
%FITDETERMINISTIC Stage 1 only: bounded nonlinear LS, with no covariance loop.
% L is the local derivative of the fitted parameters with respect to y.
values = [definition.parameters.value]';
free = find([definition.parameters.estimate]);
fit = struct('values',values,'free',free,'context',context,'definition',definition);
if isempty(free)
    [fit.mean,~] = phase_model_custom.evaluate(definition,values,coordinates,context);
    fit.J = zeros(numel(observations),0); fit.L = zeros(0,numel(observations));
    fit.assignedSensitivity=zeros(0,numel(values));
    return
end
lo = [definition.parameters(free).lower]'; hi = [definition.parameters(free).upper]';
scale = max(1,abs(values(free)));
p0 = values(free);
options = optimset('Display','off','MaxIter',1500,'MaxFunEvals',10000, ...
    'TolX',1e-8,'TolFun',1e-8,'OutputFcn', ...
    @(x,values,state) phase_model_custom.checkStop(context));
if exist('lsqnonlin','file') == 2
    [q,~,~,exitflag] = lsqnonlin(@residual,zeros(size(p0)), ...
        (lo-p0)./scale,(hi-p0)./scale,options);
else
    % Smooth bound transforms keep every trial inside the permitted domain.
    start = inverseTransform(p0,lo,hi);
    [z,~,exitflag] = fminsearch(@(z) sum(residual((transform(z,lo,hi)-p0)./scale).^2),start,options);
    q = (transform(z,lo,hi)-p0)./scale;
end
if exitflag <= 0
    error('PHASE:customFitConvergence','Custom deterministic fit did not converge. Review starting values and bounds.');
end
values(free) = p0+scale.*q;
distance = min(values(free)-lo,hi-values(free));
if any(distance < 1e-6*max(1,abs(values(free))))
    error('PHASE:customFitBoundary','An estimated parameter is on a bound; local Gaussian uncertainty is invalid. Fix that parameter or review the bounds.');
end
[fit.mean,Jall] = phase_model_custom.evaluate(definition,values,coordinates,context);
J = Jall(:,free);
if rank(J) ~= numel(free)
    error('PHASE:customIdentifiability','The custom parameters are not identifiable after temporal referencing. Remove redundant parameters.');
end
% Observed Hessian, including residual curvature for nonlinear models.
% g = J''(m-y); implicit differentiation gives dp/dy = H\J''.
H = zeros(numel(free));
for k = 1:numel(free)
    h = eps^(1/4)*max(1,abs(values(free(k))));
    distance = min(values(free(k))-lo(k),hi(k)-values(free(k)));
    if isfinite(distance), h=min(h,distance/4); end
    plus = values; minus = values;
    plus(free(k)) = plus(free(k))+h; minus(free(k)) = minus(free(k))-h;
    H(:,k) = (gradientAt(plus)-gradientAt(minus))/(2*h);
end
H = (H+H')/2;
sc = max(sqrt(sum(J.^2,1)),eps);
[~,bad] = chol(H./(sc'*sc));
if bad || rcond(H./(sc'*sc)) < 1e-11
    error('PHASE:customIdentifiability','The custom fit has an indefinite or ill-conditioned Hessian. Simplify or rescale its parameters.');
end
fit.values = values; fit.J = J; fit.L = H\J';
fit.assignedSensitivity=zeros(numel(free),numel(values));
if isfield(definition,'assignedCovariance')
    assigned=find(diag(definition.assignedCovariance)>0);
    for k=assigned(:)'
        h=eps^(1/4)*max(1,abs(values(k)));
        plus=values;minus=values;plus(k)=plus(k)+h;minus(k)=minus(k)-h;
        fit.assignedSensitivity(:,k)=-H\((gradientAt(plus)-gradientAt(minus))/(2*h));
    end
end
    function r = residual(q)
        v = values; v(free) = p0+scale.*q;
        r = phase_model_custom.evaluate(definition,v,coordinates,context)-observations;
    end
    function g = gradientAt(v)
        [m,j] = phase_model_custom.evaluate(definition,v,coordinates,context);
        g = j(:,free)'*(m-observations);
    end
end
function z = inverseTransform(p,lo,hi)
z = p;
for k=1:numel(p)
    if isfinite(lo(k)) && isfinite(hi(k)), z(k)=log((p(k)-lo(k))/(hi(k)-p(k)));
    elseif isfinite(lo(k)), z(k)=log(p(k)-lo(k));
    elseif isfinite(hi(k)), z(k)=log(hi(k)-p(k)); end
end
end
function p = transform(z,lo,hi)
p=z;
for k=1:numel(z)
    if isfinite(lo(k)) && isfinite(hi(k)), p(k)=lo(k)+(hi(k)-lo(k))/(1+exp(-z(k)));
    elseif isfinite(lo(k)), p(k)=lo(k)+exp(z(k));
    elseif isfinite(hi(k)), p(k)=hi(k)-exp(z(k)); end
end
end
