function result = modelAnomalies(displacement,uncertainty)
%MODELANOMALIES Explore deformation contrasting with the modelled AOI field.
% This is a descriptive, AOI-relative score, not a hazard probability or
% a statistical p-value. A regular 1D/2D grid samples length/area evenly.

if ~ismatrix(displacement) || isempty(displacement) || ...
        ~isequal(size(displacement),size(uncertainty))
    error('PHASE_Model_beta:invalidAnomalyInput', ...
        'Modelled displacement and uncertainty must have matching nonempty dimensions.');
end

[nPoints,nEpochs] = size(displacement);
baseline = median(displacement,1,'omitnan');
contrast = displacement-baseline;
spatialScale = 1.4826*median(abs(contrast),1,'omitnan');
uncertaintyFloor = median(uncertainty,1,'omitnan');
uncertaintyFloor(~isfinite(uncertaintyFloor) | uncertaintyFloor<0) = 0;
scale = max(spatialScale,uncertaintyFloor);
% This is only a floating-point guard, not a physical anomaly threshold.
numericFloor = 1e-9*max(1,max(abs(displacement),[],1,'omitnan'));
numericFloor(~isfinite(numericFloor)) = 1e-9;
scale = max(scale,numericFloor);
standardized = abs(contrast)./scale;
standardized(~isfinite(displacement)) = NaN;

result = struct();
result.baseline = baseline;
result.scale = scale;
result.peakDeviationMm = max(abs(contrast),[],2,'omitnan');
result.peakContrast = max(standardized,[],2,'omitnan');
result.persistence = sum(standardized>=3,2)./max(1,sum(isfinite(standardized),2));
result.class = repmat({'typical'},nPoints,1);

% A short transient needs two consecutive modelled epochs above threshold.
for point = 1:nPoints
    above = isfinite(standardized(point,:)) & standardized(point,:)>=3;
    if nEpochs>=2 && any(above(1:end-1) & above(2:end))
        if result.persistence(point)>=0.2
            result.class{point} = 'persistent';
        else
            result.class{point} = 'transient';
        end
    end
end
end
