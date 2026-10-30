function K=kernel(a,b,settings)
%KERNEL Original separable analytical covariance fitted to empirical products.
if strcmp(settings.dimension,'1D'), ds=abs(a(:,3)-b(:,3)');
else, ds=hypot(a(:,1)-b(:,1)',a(:,2)-b(:,2)'); end
dt=abs(a(:,4)-b(:,4)');
K=settings.timeFunction(settings.timeParameters,dt).* ...
    settings.spaceFunction(settings.spaceParameters,ds)/settings.timeParameters(1);
end
