function ui = configToUi(cfg)
%CONFIGTOUI Convert MATLAB values to JSON-safe form controls.

ui = struct();
ui.constellation = ternary(strcmp(cfg.constellation, 'CSK'), 'COSMO-SkyMed', 'Sentinel-1');
textFields = {'python','master_date','polarisation','dem_name','dem_file', ...
    'dem_name_coreg','dem_file_coreg','dem_resampling','gptbin_path','cache'};
for k = 1:numel(textFields)
    name = textFields{k};
    ui.(name) = char(string(cfg.(name)));
end
numericFields = {'lon_min','lat_min','lon_max','lat_max','first_step', ...
    'epsg_code','cpu','num_gcp'};
for k = 1:numel(numericFields)
    name = numericFields{k};
    ui.(name) = num2str(cfg.(name), 15);
end
ui.update_processed_data = logical(cfg.update_processed_data);
ui.auto_master = logical(cfg.auto_master);
ui.process_master = logical(cfg.process_master);
ui.remove_slaves_after_processing = logical(cfg.remove_slaves_after_processing);
ui.remove_split_after_processing = logical(cfg.remove_split_after_processing);
ui.remove_coreg_after_processing = logical(cfg.remove_coreg_after_processing);
ui.remove_ifg_after_processing = logical(cfg.remove_ifg_after_processing);
ui.auto_epsg = logical(cfg.auto_epsg);
ui.generate_coherence = logical(cfg.generate_coherence);
ui.generate_lia = logical(cfg.generate_lia);
end

function out = ternary(condition, yesValue, noValue)
if condition, out = yesValue; else, out = noValue; end
end
