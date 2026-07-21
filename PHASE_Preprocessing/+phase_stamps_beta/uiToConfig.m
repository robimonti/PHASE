function cfg = uiToConfig(ui, base)
%UITOCONFIG Normalize values received from the HTML form.

if nargin < 2 || isempty(base)
    base = phase_stamps_beta.defaultConfig();
end
cfg = base;

cfg.stamps_preparation = double(~logicalValue(ui, 'prepare_data'));
cfg.train_flag = double(~logicalValue(ui, 'train_enabled'));

textFields = {'installation_folder','project_path','master_date','export_name', ...
    'utc_time','filter_weighting','quick_est_gamma_flag','small_baseline_flag', ...
    'select_method','weed_neighbours','weed_zero_elevation','unwrap_method', ...
    'unwrap_spatial_cost_func_flag','unwrap_prefilter_flag','unwrap_patch_phase', ...
    'unwrap_la_error_flag','unwrap_hold_good_values','subtr_tropo','tropo_method', ...
    'select_reest_gamma_flag','scla_deramp','scla_method','scn_kriging_flag', ...
    'ph_output','stamps_first_step','stamps_last_step'};
for k = 1:numel(textFields)
    name = textFields{k};
    if isfield(ui, name)
        cfg.(name) = char(string(ui.(name)));
    end
end
cfg.master_date = regexprep(cfg.master_date, '[^0-9]', '');

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
    if isfield(ui, name)
        cfg.(name) = numericValue(ui.(name), name);
    end
end

if isfield(ui, 'initial_date')
    try
        d = datetime(char(string(ui.initial_date)), 'InputFormat', 'yyyy-MM-dd');
        cfg.year_0 = year(d);
        cfg.month_0 = month(d);
        cfg.day_0 = day(d);
    catch
        error('PHASE_StaMPS_beta:invalidDate', ...
            'Initial acquisition date must use YYYY-MM-DD.');
    end
end

cfg.drop_ifg_index = vectorValue(ui, 'drop_ifg_index');
cfg.scla_drop_index = vectorValue(ui, 'scla_drop_index');
cfg.ref_centre_lonlat = [numericValue(ui.ref_centre_lon, 'ref_centre_lon'), ...
    numericValue(ui.ref_centre_lat, 'ref_centre_lat')];
cfg.ref_centre_lonlat_w = [numericValue(ui.ref_centre_lon_w, 'ref_centre_lon_w'), ...
    numericValue(ui.ref_centre_lat_w, 'ref_centre_lat_w')];
end

function value = logicalValue(ui, name)
if ~isfield(ui, name)
    value = false;
    return
end
raw = ui.(name);
if islogical(raw)
    value = raw;
elseif isnumeric(raw)
    value = raw ~= 0;
else
    value = any(strcmpi(char(string(raw)), {'true','1','yes','y','on'}));
end
end

function value = numericValue(raw, name)
if isnumeric(raw) && isscalar(raw)
    value = double(raw);
else
    value = str2double(strtrim(char(string(raw))));
end
if isempty(value) || ~isscalar(value) || isnan(value)
    error('PHASE_StaMPS_beta:invalidNumber', ...
        '%s must be a valid numeric value.', name);
end
end

function out = vectorValue(ui, name)
if ~isfield(ui, name)
    out = '[]';
    return
end
raw = strtrim(char(string(ui.(name))));
if isempty(raw) || strcmp(raw, '[]')
    out = '[]';
    return
end
raw = regexprep(raw, '[\[\],;]', ' ');
tokens = strsplit(strtrim(raw));
values = cellfun(@str2double, tokens);
if any(isnan(values)) || any(values < 1) || any(mod(values, 1) ~= 0)
    error('PHASE_StaMPS_beta:invalidIndexVector', ...
        '%s must be [] or a list of positive integer indices.', name);
end
out = strtrim(sprintf('%d ', values));
end
