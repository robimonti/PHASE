function validateOptions(config)
%VALIDATEOPTIONS Validate custom mean and uncertainty controls.
if ~ismember(config.customResidualMean,{'model','constant'})
    error('PHASE:customMean','Select model-only diagnostics or an estimated constant residual offset.');
end
if ~ismember(config.customCoordinates,{'project','lonlat'})
    error('PHASE:customCoordinates','Select project or longitude/latitude coordinates.');
end
for name={'customOriginX','customOriginY'}
    value=config.(name{1});
    if ~isnumeric(value)||~isscalar(value)||(~isnan(value)&&~isfinite(value))
        error('PHASE:customOrigin','Origins must be finite or empty for automatic centering.');
    end
end
if ~isscalar(config.customTimeOffset)||~isfinite(config.customTimeOffset)
    error('PHASE:customTime','The model time offset must be finite.');
end
limit=config.customMaxObservations;
if ~isnumeric(limit)||~isscalar(limit)||~isfinite(limit)||limit<5||limit~=round(limit)
    error('PHASE:customOption','The dense observation limit must be an integer of at least 5.');
end
end
