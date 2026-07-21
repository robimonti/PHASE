function ui = configToUi(cfg)
%CONFIGTOUI Convert MAT-compatible values to JSON-safe form controls.

ui = struct();
ui.prepare_data = logical(cfg.stamps_preparation == 0);
ui.train_enabled = logical(cfg.train_flag == 0);

textFields = {'installation_folder','project_path','master_date','export_name', ...
    'utc_time','filter_weighting','quick_est_gamma_flag','small_baseline_flag', ...
    'select_method','weed_neighbours','weed_zero_elevation','unwrap_method', ...
    'unwrap_spatial_cost_func_flag','unwrap_prefilter_flag','unwrap_patch_phase', ...
    'unwrap_la_error_flag','unwrap_hold_good_values','subtr_tropo','tropo_method', ...
    'select_reest_gamma_flag','scla_deramp','scla_method','scn_kriging_flag', ...
    'ph_output','stamps_first_step','stamps_last_step'};
for k = 1:numel(textFields)
    name = textFields{k};
    ui.(name) = char(string(cfg.(name)));
end

numericFields = {'amplitude_threshold','time_span','weed_time_win', ...
    'unwrap_time_win','scn_time_win','n_cores','heading','lambda','max_topo_err', ...
    'filter_grid_size','gamma_max_iterations','gamma_change_convergence', ...
    'gamma_stdev_reject','clap_win','clap_alpha','clap_beta', ...
    'clap_low_pass_wavelength','percent_rand','density_rand','weed_standard_dev', ...
    'weed_max_noise','merge_resample_size','merge_standard_dev','unwrap_grid_size', ...
    'unwrap_gold_n_win','unwrap_gold_alpha','unwrap_alpha','scn_wavelength', ...
    'ref_radius','ref_velocity','plot_s','ref_radius_w'};
for k = 1:numel(numericFields)
    name = numericFields{k};
    ui.(name) = scalarText(cfg.(name));
end

ui.initial_date = sprintf('%04d-%02d-%02d', cfg.year_0, cfg.month_0, cfg.day_0);
ui.drop_ifg_index = vectorText(cfg.drop_ifg_index);
ui.scla_drop_index = vectorText(cfg.scla_drop_index);
ui.ref_centre_lon = scalarText(cfg.ref_centre_lonlat(1));
ui.ref_centre_lat = scalarText(cfg.ref_centre_lonlat(2));
ui.ref_centre_lon_w = scalarText(cfg.ref_centre_lonlat_w(1));
ui.ref_centre_lat_w = scalarText(cfg.ref_centre_lonlat_w(2));
end

function out = scalarText(value)
if ischar(value) || isstring(value)
    out = char(string(value));
elseif isinf(value)
    if value > 0, out = 'Inf'; else, out = '-Inf'; end
elseif isnan(value)
    out = 'NaN';
else
    out = num2str(value, 15);
end
end

function out = vectorText(value)
if ischar(value) || isstring(value)
    out = char(string(value));
elseif isempty(value)
    out = '[]';
else
    out = strtrim(regexprep(mat2str(value), '[\[\],;]', ' '));
end
end
