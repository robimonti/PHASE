function cfg = uiToConfig(ui, base)
%UITOCONFIG Normalize values received from the HTML form.

if nargin < 2 || isempty(base), base = phase_preprocessing_beta.defaultConfig(); end
cfg = base;
if isfield(ui, 'constellation')
    value = upper(char(string(ui.constellation)));
    cfg.constellation = ternary(startsWith(value, 'COSMO') || strcmp(value, 'CSK'), 'CSK', 'SEN');
end

textFields = {'python','master_date','polarisation','dem_name','dem_file', ...
    'dem_name_coreg','dem_file_coreg','dem_resampling','gptbin_path','cache'};
for k = 1:numel(textFields)
    name = textFields{k};
    if isfield(ui, name), cfg.(name) = char(string(ui.(name))); end
end
cfg.master_date = regexprep(cfg.master_date, '[^0-9]', '');

numericFields = {'lon_min','lat_min','lon_max','lat_max','first_step', ...
    'epsg_code','cpu','num_gcp'};
for k = 1:numel(numericFields)
    name = numericFields{k};
    if isfield(ui, name), cfg.(name) = numericValue(ui.(name), name); end
end
logicalFields = {'update_processed_data','auto_master','process_master', ...
    'remove_slaves_after_processing','remove_split_after_processing', ...
    'remove_coreg_after_processing','remove_ifg_after_processing', ...
    'auto_epsg','generate_coherence','generate_lia'};
for k = 1:numel(logicalFields)
    name = logicalFields{k};
    if isfield(ui, name), cfg.(name) = logicalValue(ui.(name)); end
end
if cfg.first_step > 1
    % Resuming the slave pipeline requires the existing master products.
    % Reprocessing the master would reset those dependencies and is therefore
    % both unnecessary and unsafe.
    cfg.process_master = false;
end
if cfg.auto_epsg
    cfg.epsg_code = phase_preprocessing_beta.estimateEpsg( ...
        cfg.lon_min,cfg.lat_min,cfg.lon_max,cfg.lat_max);
end
end

function value = numericValue(raw, name)
if isnumeric(raw) && isscalar(raw), value = double(raw); else, value = str2double(strtrim(char(string(raw)))); end
if isempty(value) || ~isscalar(value) || ~isfinite(value)
    error('PHASE_Preprocessing_beta:invalidNumber', '%s must be a valid number.', name);
end
end

function value = logicalValue(raw)
if islogical(raw), value = raw;
elseif isnumeric(raw), value = raw ~= 0;
else, value = any(strcmpi(char(string(raw)), {'true','1','yes','y','on'}));
end
end

function out = ternary(condition, yesValue, noValue)
if condition, out = yesValue; else, out = noValue; end
end
