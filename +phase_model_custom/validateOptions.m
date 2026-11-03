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
if isfield(config,'customResidualField') && ~ismember(config.customResidualField,{'automatic','common'})
    error('PHASE:customResidualField','Unknown residual deterministic field policy.');
end
limit=config.customMemoryGB;
if ~isnumeric(limit)||~isscalar(limit)||~isfinite(limit)||limit<0.1
    error('PHASE:customOption','The STC memory budget must be at least 0.1 GiB.');
end
end
