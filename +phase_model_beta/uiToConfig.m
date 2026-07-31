function config = uiToConfig(payload)
%UITOCONFIG Convert HTML form values into the typed Model configuration.

config = phase_model_beta.defaultConfig();
definition = phase_model_beta.schema();
for k = 1:numel(definition.fields)
    meta = definition.fields(k);
    id = meta.id;
    if ~isfield(payload,id), continue; end
    raw = payload.(id);
    switch meta.type
        case 'boolean'
            config.(id) = logical(raw);
        case 'number'
            if isempty(raw)
                config.(id) = NaN;
            else
                value = double(raw);
                if ~isscalar(value) || ~isfinite(value)
                    error('PHASE_Model_beta:invalidNumber', ...
                        '%s must be a finite scalar or empty.',meta.label);
                end
                config.(id) = value;
            end
        case 'date'
            try
                config.(id) = datetime(char(string(raw)), ...
                    'InputFormat','yyyy-MM-dd');
            catch
                error('PHASE_Model_beta:invalidDate', ...
                    '%s must be a valid calendar date.',meta.label);
            end
        otherwise
            config.(id) = char(string(raw));
    end
end

if isfield(payload,'aoi_polygon_lonlat')
    config.aoi_polygon_lonlat = numericMatrix(payload.aoi_polygon_lonlat);
end

covariancePairs = {
    'dtCov_method_STC1D','dtCov_STC1D'
    'dsCov_method_STC1D','dsCov_STC1D'
    'dtCov_method_STC2D','dtCov_STC2D'
    'dsCov_method_STC2D','dsCov_STC2D'
};
for k = 1:size(covariancePairs,1)
    if strcmp(config.(covariancePairs{k,1}),'auto')
        config.(covariancePairs{k,2}) = NaN;
    end
end
end

function value = numericMatrix(raw)
if isempty(raw)
    value = zeros(0,2);
elseif iscell(raw)
    try
        if isvector(raw) && all(cellfun(@(row) isnumeric(row) && numel(row)==2,raw))
            value = vertcat(raw{:});
        else
            value = cell2mat(raw);
        end
    catch
        error('PHASE_Model_beta:invalidPolygon','The AOI polygon coordinates are invalid.');
    end
else
    value = double(raw);
end
value = squeeze(value);
if size(value,2) ~= 2 && size(value,1) == 2, value = value.'; end
if ~isempty(value) && (size(value,2) ~= 2 || any(~isfinite(value(:))))
    error('PHASE_Model_beta:invalidPolygon', ...
        'The AOI polygon must contain finite [longitude, latitude] rows.');
end
end
