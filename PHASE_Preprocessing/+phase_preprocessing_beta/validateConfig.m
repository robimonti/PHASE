function [errors, warnings] = validateConfig(cfg, forStart)
%VALIDATECONFIG Guard the processing inputs before saving or execution.

if nargin < 2, forStart = false; end
errors = {}; warnings = {};
if ~any(strcmp(cfg.constellation, {'SEN','CSK'})), errors{end+1} = 'Select Sentinel-1 or COSMO-SkyMed.'; end
if isempty(strtrim(cfg.python)), errors{end+1} = 'Python executable/environment cannot be empty.'; end
if isempty(strtrim(cfg.gptbin_path)), errors{end+1} = 'SNAP gpt path cannot be empty.'; end
if cfg.lon_min >= cfg.lon_max, errors{end+1} = 'Minimum longitude must be smaller than maximum longitude.'; end
if cfg.lat_min >= cfg.lat_max, errors{end+1} = 'Minimum latitude must be smaller than maximum latitude.'; end
if cfg.lon_min < -180 || cfg.lon_max > 180, errors{end+1} = 'Longitude must stay between -180 and 180 degrees.'; end
if cfg.lat_min < -90 || cfg.lat_max > 90, errors{end+1} = 'Latitude must stay between -90 and 90 degrees.'; end
if ~cfg.auto_master && isempty(regexp(cfg.master_date, '^\d{8}$', 'once'))
    errors{end+1} = 'Manual master date must use YYYYMMDD.';
end
if cfg.first_step < 1 || cfg.first_step > 6 || mod(cfg.first_step, 1) ~= 0
    errors{end+1} = 'First preprocessing step must be an integer from 1 to 6.';
end
if cfg.first_step == 6 && ~(cfg.generate_coherence || cfg.generate_lia)
    errors{end+1} = 'Step 6 requires coherence, LIA, or both to be enabled.';
end
if cfg.first_step > 1
    warnings{end+1} = sprintf([ ...
        'Resuming from slave step %d: master processing is disabled and ', ...
        'the existing master products will be reused.'],cfg.first_step);
end
if cfg.cpu < 1 || mod(cfg.cpu, 1) ~= 0, errors{end+1} = 'CPU cores must be a positive integer.'; end
if cfg.epsg_code < 1 || mod(cfg.epsg_code, 1) ~= 0, errors{end+1} = 'EPSG code must be a positive integer.'; end
if isempty(errors) || ~any(contains(errors,'EPSG'))
    try
        projcrs(cfg.epsg_code);
    catch
        errors{end+1} = sprintf('EPSG:%d is not a supported projected CRS.',cfg.epsg_code);
    end
end
if strcmp(cfg.constellation, 'CSK') && (cfg.num_gcp < 1 || mod(cfg.num_gcp, 1) ~= 0)
    errors{end+1} = 'COSMO-SkyMed GCP count must be a positive integer.';
end
if strcmp(cfg.dem_name, 'External DEM') && isempty(strtrim(cfg.dem_file))
    errors{end+1} = 'Select the external interferogram DEM file.';
end
if strcmp(cfg.constellation, 'SEN') && strcmp(cfg.dem_name_coreg, 'External DEM') && isempty(strtrim(cfg.dem_file_coreg))
    errors{end+1} = 'Select the external coregistration DEM file.';
end
if forStart
    try
        [~,gptInfo] = phase_preprocessing_beta.resolveGpt(cfg.gptbin_path);
        if ~gptInfo.exists
            warnings{end+1} = ['SNAP gpt executable was not found on this computer: ' cfg.gptbin_path];
        end
    catch ME
        warnings{end+1} = ['SNAP gpt executable is not usable: ' ME.message];
    end
end
end
