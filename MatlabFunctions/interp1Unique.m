function values = interp1Unique(samplePoints,sampleValues,queryPoints,varargin)
%INTERP1UNIQUE One-dimensional interpolation after consolidating duplicates.
%
% interp1 delegates to griddedInterpolant, which requires unique sample
% coordinates. Repeated acquisition times or empirical-covariance lags are
% consolidated using their mean value before preserving the requested
% interpolation method and extrapolation behaviour.

originalQuerySize = size(queryPoints);
x = double(samplePoints(:));
y = double(sampleValues(:));
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if isempty(x)
    values = NaN(originalQuerySize);
    return
end

[x,order] = sort(x);
y = y(order);
[uniqueX,~,groups] = unique(x);
if numel(uniqueX) < numel(x)
    uniqueY = accumarray(groups,y,[],@(items) mean(items,'omitnan'));
else
    uniqueY = y;
end

if numel(uniqueX) == 1
    values = repmat(uniqueY,originalQuerySize);
    return
end

values = interp1(uniqueX,uniqueY,queryPoints,varargin{:});
end
