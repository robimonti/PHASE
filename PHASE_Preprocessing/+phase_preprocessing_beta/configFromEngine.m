function cfg = configFromEngine(engine, base)
%CONFIGFROMENGINE Import edits made in the advanced MATLAB workspace.

cfg = base;
cfg.constellation = char(string(engine.constellation));
if strcmp(cfg.constellation, 'SEN')
    cfg.python = char(string(engine.python_SEN));
    cfg.update_processed_data = logical(engine.update_processed_data_SEN);
    cfg.master_date = dateText(engine.master_date_SEN);
    cfg.auto_master = logical(engine.auto_master_SEN);
    cfg.process_master = engine.master_processing_SEN == 0;
    cfg.polarisation = char(string(engine.polarisation_SEN));
    cfg.lon_min = engine.lon_min_SEN; cfg.lon_max = engine.lon_max_SEN;
    cfg.lat_min = engine.lat_min_SEN; cfg.lat_max = engine.lat_max_SEN;
    cfg.remove_slaves_after_processing = engine.slaves_removal_SEN == 0;
    cfg.dem_name = char(string(engine.dem_name_SEN));
    cfg.dem_file = char(string(engine.dem_file_SEN));
    cfg.dem_name_coreg = char(string(engine.dem_name_coreg_SEN));
    cfg.dem_file_coreg = char(string(engine.dem_file_coreg_SEN));
    cfg.dem_resampling = char(string(engine.dem_resampling_SEN));
    cfg.first_step = numeric(engine.first_step_SEN, cfg.first_step);
    cfg.generate_coherence = engine.coherence_tc_SEN == 0;
    cfg.generate_lia = cfg.generate_coherence;
    cfg.epsg_code = engine.epsg_code_SEN;
    cfg.gptbin_path = char(string(engine.gptbin_path_SEN));
    cfg.cpu = engine.cpu_SEN; cfg.cache = char(string(engine.cache_SEN));
else
    cfg.python = char(string(engine.python_CSK));
    cfg.master_date = dateText(engine.master_date_CSK);
    cfg.auto_master = logical(engine.auto_master_CSK);
    cfg.process_master = engine.master_processing_CSK == 0;
    cfg.lon_min = engine.lon_min_CSK; cfg.lon_max = engine.lon_max_CSK;
    cfg.lat_min = engine.lat_min_CSK; cfg.lat_max = engine.lat_max_CSK;
    cfg.remove_slaves_after_processing = engine.slaves_removal_CSK == 0;
    cfg.dem_name = char(string(engine.dem_name_CSK));
    cfg.dem_file = char(string(engine.dem_file_CSK));
    cfg.num_gcp = engine.num_gcp_CSK;
    cfg.first_step = numeric(engine.first_step_CSK, cfg.first_step);
    cfg.generate_coherence = engine.coherence_tc_CSK == 0;
    cfg.generate_lia = cfg.generate_coherence;
    cfg.epsg_code = engine.epsg_code_CSK;
    cfg.gptbin_path = char(string(engine.gptbin_path_CSK));
    cfg.cpu = engine.cpu_CSK; cfg.cache = char(string(engine.cache_CSK));
end
cfg.master_date = regexprep(cfg.master_date, '[^0-9]', '');
end

function out = dateText(value)
if isdatetime(value)
    out = char(datetime(value, 'Format', 'yyyyMMdd'));
else
    out = char(string(value));
end
end

function out = numeric(value, fallback)
if isnumeric(value), out = double(value); else, out = str2double(char(string(value))); end
if ~isscalar(out) || ~isfinite(out), out = fallback; end
end
