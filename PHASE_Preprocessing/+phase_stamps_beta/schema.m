function out = schema()
%SCHEMA UI metadata kept outside both the HTML and the processing engine.

items = repmat(item('', '', '', '', '', '', {}, false), 0, 1);

items(end+1) = item('prepare_data', 'Prepare data for StaMPS', 'project', 'toggle', '', ...
    'Run mt_prep_snap before StaMPS processing.', {}, false);
items(end+1) = item('installation_folder', 'StaMPS installation folder', 'project', 'path', '', ...
    'Root containing matlab, matlab_compat and the native binaries.', {}, false);
items(end+1) = item('project_path', 'PHASE project folder', 'project', 'path', '', ...
    'PHASE root containing PHASE_Preprocessing/INSAR_<master date>.', {}, false);
items(end+1) = item('amplitude_threshold', 'Amplitude threshold', 'project', 'number', '', ...
    'Threshold passed to mt_prep_snap.', {}, false);
items(end+1) = item('master_date', 'Master date', 'project', 'text', 'YYYYMMDD', ...
    'Reference acquisition date.', {}, false);
items(end+1) = item('export_name', 'Export filename', 'project', 'text', '', ...
    'Base name used for XLSX/CSV exports.', {}, false);

items(end+1) = item('initial_date', 'Initial acquisition date', 'global', 'date', '', ...
    'First acquisition in the stack.', {}, false);
items(end+1) = item('time_span', 'Dataset time span', 'global', 'number', 'days', ...
    'Automatically calculated from first to last acquisition.', {}, false);
items(end+1) = item('utc_time', 'Satellite UTC time', 'global', 'text', 'HH:MM', ...
    'Rounded center time used by TRAIN/GACOS.', {}, false);
items(end+1) = item('heading', 'Heading', 'global', 'number', 'deg', ...
    'Satellite heading extracted from the SNAP parameter file.', {}, false);
items(end+1) = item('lambda', 'Radar wavelength', 'global', 'number', 'm', ...
    'Computed from radar frequency.', {}, false);
items(end+1) = item('n_cores', 'Number of cores', 'global', 'number', '', ...
    'Requested parallel workers.', {}, false);

items(end+1) = item('train_enabled', 'TRAIN atmospheric correction', 'step1', 'toggle', '', ...
    'Enable TRAIN and the selected tropospheric model.', {}, false);
items(end+1) = item('stamps_first_step', 'StaMPS first step', 'step1', 'choice', '', ...
    'First processing step to execute.', {'1','2','3','4','5','6'}, false);
items(end+1) = item('stamps_last_step', 'StaMPS last step', 'step1', 'choice', '', ...
    'Final processing step to execute.', {'7','8'}, false);
items(end+1) = item('max_topo_err', 'Maximum topographic error', 'step1', 'number', 'm', ...
    'Maximum DEM error used in phase analysis.', {}, false);

items(end+1) = item('filter_grid_size', 'Filter grid size', 'step2', 'number', 'm', '', {}, false);
items(end+1) = item('filter_weighting', 'Filter weighting', 'step2', 'text', '', '', {}, true);
items(end+1) = item('gamma_max_iterations', 'Gamma maximum iterations', 'step2', 'number', '', '', {}, false);
items(end+1) = item('gamma_change_convergence', 'Gamma convergence change', 'step2', 'number', '', '', {}, true);
items(end+1) = item('gamma_stdev_reject', 'Gamma stdev reject', 'step2', 'number', '', '', {}, true);
items(end+1) = item('quick_est_gamma_flag', 'Quick gamma estimate', 'step2', 'choice', '', '', {'y','n'}, false);
items(end+1) = item('small_baseline_flag', 'Small baseline mode', 'step2', 'choice', '', '', {'n','y'}, false);
items(end+1) = item('clap_win', 'CLAP window', 'step2', 'number', '', '', {}, false);
items(end+1) = item('clap_alpha', 'CLAP alpha', 'step2', 'number', '', '', {}, true);
items(end+1) = item('clap_beta', 'CLAP beta', 'step2', 'number', '', '', {}, true);
items(end+1) = item('clap_low_pass_wavelength', 'CLAP low-pass wavelength', 'step2', 'number', 'm', '', {}, true);

items(end+1) = item('select_method', 'Selection method', 'step3', 'choice', '', '', {'PERCENT','DENSITY'}, false);
items(end+1) = item('percent_rand', 'Random percentage', 'step3', 'number', '%', ...
    'Used when selection method is PERCENT.', {}, false);
items(end+1) = item('density_rand', 'Random density', 'step3', 'number', '', ...
    'Used when selection method is DENSITY.', {}, false);

items(end+1) = item('weed_time_win', 'Weeding temporal window', 'step4', 'number', 'days', '', {}, false);
items(end+1) = item('weed_standard_dev', 'Weeding standard deviation', 'step4', 'number', '', '', {}, false);
items(end+1) = item('weed_neighbours', 'Weed neighbours', 'step4', 'choice', '', '', {'y','n'}, false);
items(end+1) = item('weed_zero_elevation', 'Weed zero elevation', 'step4', 'choice', '', '', {'n','y'}, true);
items(end+1) = item('weed_max_noise', 'Maximum noise', 'step4', 'number', '', 'Inf disables this threshold.', {}, true);

items(end+1) = item('merge_resample_size', 'Merge resample size', 'step5', 'number', '', '', {}, false);
items(end+1) = item('merge_standard_dev', 'Merge standard deviation', 'step5', 'number', '', 'Inf disables this threshold.', {}, true);

items(end+1) = item('unwrap_grid_size', 'Unwrap grid size', 'step6', 'number', 'm', '', {}, false);
items(end+1) = item('unwrap_gold_n_win', 'Goldstein window', 'step6', 'number', '', '', {}, false);
items(end+1) = item('unwrap_time_win', 'Unwrap temporal window', 'step6', 'number', 'days', '', {}, false);
items(end+1) = item('unwrap_method', 'Unwrap method', 'step6', 'text', '', '', {}, false);
items(end+1) = item('unwrap_gold_alpha', 'Goldstein alpha', 'step6', 'number', '', '', {}, true);
items(end+1) = item('unwrap_alpha', 'Unwrap alpha', 'step6', 'number', '', '', {}, true);
items(end+1) = item('unwrap_spatial_cost_func_flag', 'Spatial cost function', 'step6', 'choice', '', '', {'n','y'}, true);
items(end+1) = item('unwrap_prefilter_flag', 'Prefilter', 'step6', 'choice', '', '', {'y','n'}, true);
items(end+1) = item('unwrap_patch_phase', 'Patch phase', 'step6', 'choice', '', '', {'n','y'}, true);
items(end+1) = item('unwrap_la_error_flag', 'Look-angle error', 'step6', 'choice', '', '', {'y','n'}, true);
items(end+1) = item('unwrap_hold_good_values', 'Hold good values', 'step6', 'choice', '', '', {'y','n'}, true);

items(end+1) = item('subtr_tropo', 'Subtract troposphere', 'step7', 'choice', '', '', {'y','n'}, false);
items(end+1) = item('tropo_method', 'Tropospheric method', 'step7', 'choice', '', '', {'a_gacos','a_linear'}, false);
items(end+1) = item('select_reest_gamma_flag', 'Re-estimate gamma', 'step7', 'choice', '', '', {'y','n'}, true);
items(end+1) = item('drop_ifg_index', 'Drop IFG indices', 'step7', 'vector', '', '[] or a space-separated list such as 1 3.', {}, true);
items(end+1) = item('scla_deramp', 'SCLA deramp', 'step7', 'choice', '', '', {'y','n'}, false);
items(end+1) = item('scla_method', 'SCLA method', 'step7', 'text', '', '', {}, false);
items(end+1) = item('scla_drop_index', 'SCLA drop indices', 'step7', 'vector', '', '[] or a space-separated list.', {}, true);
items(end+1) = item('scn_wavelength', 'Atmospheric spatial wavelength', 'step8', 'number', 'm', '', {}, false);
items(end+1) = item('scn_time_win', 'Atmospheric temporal window', 'step8', 'number', 'days', '', {}, false);
items(end+1) = item('scn_kriging_flag', 'SCN kriging', 'step8', 'choice', '', '', {'n','y'}, true);

items(end+1) = item('ph_output', 'Phase output', 'export', 'choice', '', '', {'unwrapped','wrapped'}, false);
items(end+1) = item('ref_centre_lon', 'Reference longitude', 'export', 'number', 'deg', '', {}, false);
items(end+1) = item('ref_centre_lat', 'Reference latitude', 'export', 'number', 'deg', '', {}, false);
items(end+1) = item('ref_radius', 'Reference radius', 'export', 'number', 'm', '', {}, false);
items(end+1) = item('ref_velocity', 'Reference velocity', 'export', 'number', 'mm/year', '', {}, false);
items(end+1) = item('plot_s', 'Plot marker size', 'export', 'number', '', '', {}, true);
items(end+1) = item('ref_centre_lon_w', 'Wrapped reference longitude', 'export', 'number', 'deg', '', {}, true);
items(end+1) = item('ref_centre_lat_w', 'Wrapped reference latitude', 'export', 'number', 'deg', '', {}, true);
items(end+1) = item('ref_radius_w', 'Wrapped reference radius', 'export', 'number', 'm', '', {}, true);

groups = struct( ...
    'id', {'project','global','step1','step2','step3','step4','step5','step6','step7','step8','export'}, ...
    'title', {'Project','Global variables','Step 1','Step 2','Step 3','Step 4','Step 5','Step 6','Step 7','Step 8','Export'}, ...
    'subtitle', { ...
        'Paths, preparation and output identity', ...
        'Acquisition geometry and dataset metadata', ...
        'Load data and processing range', ...
        'Estimate phase noise', ...
        'PS selection', ...
        'PS weeding', ...
        'Merge patches', ...
        'Phase unwrapping', ...
        'Spatially correlated look-angle error', ...
        'Atmospheric filtering', ...
        'Reference area and output products'});

out = struct('groups', groups, 'items', items);
end

function s = item(id, label, group, type, unit, help, options, advanced)
s = struct('id', id, 'label', label, 'group', group, 'type', type, ...
    'unit', unit, 'help', help, 'options', {options}, 'advanced', advanced);
end
