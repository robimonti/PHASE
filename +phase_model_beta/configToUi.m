function value = configToUi(config)
%CONFIGTOUI Convert MATLAB-only values into JSON-friendly values.

value = config;
value.t0IN = char(datetime(config.t0IN,'Format','yyyy-MM-dd'));
names = fieldnames(value);
for k = 1:numel(names)
    item = value.(names{k});
    if isnumeric(item) && isscalar(item) && isnan(item)
        value.(names{k}) = [];
    end
end
end
