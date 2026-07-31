function [cfg, info] = loadConfig(rootDir)
%LOADCONFIG Read input_preprocessing.mat without relying on App Designer.

cfg = phase_preprocessing_beta.defaultConfig();
pathValue = fullfile(rootDir, 'PHASE_Preprocessing', 'input_preprocessing.mat');
info = struct('path', pathValue, 'exists', false, 'loadedFields', {{}});
if exist(pathValue, 'file') ~= 2
    return
end

data = load(pathValue);
info.exists = true;
direct = {'python','master_date','polarisation','lon_min','lat_min','lon_max', ...
    'lat_max','dem_name','dem_file','dem_name_coreg','dem_file_coreg', ...
    'dem_resampling','epsg_code','gptbin_path','cpu','cache','num_gcp'};
for k = 1:numel(direct)
    name = direct{k};
    if isfield(data, name)
        cfg.(name) = normalizeScalar(data.(name));
        info.loadedFields{end+1} = name; %#ok<AGROW>
    end
end

if isfield(data, 'constellation')
    value = upper(char(string(data.constellation)));
    if startsWith(value, 'CSK') || contains(value, 'COSMO')
        cfg.constellation = 'CSK';
    else
        cfg.constellation = 'SEN';
    end
end
if isfield(data, 'auto_master'), cfg.auto_master = logical(data.auto_master); end
if isfield(data, 'master_processing'), cfg.process_master = data.master_processing == 0; end
if isfield(data, 'slaves_removal'), cfg.remove_slaves_after_processing = data.slaves_removal == 0; end
if isfield(data, 'remove_split_after_processing')
    cfg.remove_split_after_processing = logical(data.remove_split_after_processing);
end
if isfield(data, 'remove_coreg_after_processing')
    cfg.remove_coreg_after_processing = logical(data.remove_coreg_after_processing);
end
if isfield(data, 'remove_ifg_after_processing')
    cfg.remove_ifg_after_processing = logical(data.remove_ifg_after_processing);
end
if isfield(data, 'auto_epsg')
    cfg.auto_epsg = logical(data.auto_epsg);
else
    % Existing PHASE configurations used an explicitly chosen EPSG code.
    cfg.auto_epsg = false;
end
legacyProducts = isfield(data, 'coherence_tc') && data.coherence_tc == 0;
if isfield(data, 'generate_coherence')
    cfg.generate_coherence = logical(data.generate_coherence);
else
    cfg.generate_coherence = legacyProducts;
end
if isfield(data, 'generate_lia')
    cfg.generate_lia = logical(data.generate_lia);
else
    cfg.generate_lia = legacyProducts;
end
if isfield(data, 'first_step')
    value = str2double(char(string(data.first_step)));
    if isfinite(value), cfg.first_step = value; end
end
if cfg.first_step > 1
    cfg.process_master = false;
end
cfg.master_date = regexprep(char(string(cfg.master_date)), '[^0-9]', '');
if cfg.auto_epsg
    cfg.epsg_code = phase_preprocessing_beta.estimateEpsg( ...
        cfg.lon_min,cfg.lat_min,cfg.lon_max,cfg.lat_max);
end
end

function value = normalizeScalar(value)
if isdatetime(value) && isscalar(value)
    value = char(datetime(value, 'Format', 'yyyyMMdd'));
end
if isstring(value) && isscalar(value), value = char(value); end
end
