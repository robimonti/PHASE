function [errors, warnings] = validateConfig(cfg, forStart)
%VALIDATECONFIG Validate values before Save or Start.

if nargin < 2, forStart = false; end
errors = {};
warnings = {};

if numel(cfg.master_date) ~= 8 || any(~isstrprop(cfg.master_date, 'digit'))
    errors{end+1} = 'Master date must contain exactly eight digits (YYYYMMDD).';
end
if isempty(regexp(cfg.utc_time, '^([01]\d|2[0-3]):[0-5]\d$', 'once'))
    errors{end+1} = 'UTC time must use HH:MM (00:00 to 23:59).';
end
if ~any(strcmp(cfg.stamps_first_step, {'1','2','3','4','5','6'}))
    errors{end+1} = 'StaMPS first step must be between 1 and 6.';
end
if ~any(strcmp(cfg.stamps_last_step, {'7','8'}))
    errors{end+1} = 'StaMPS last step must be 7 or 8.';
end
if ~any(strcmp(cfg.ph_output, {'wrapped','unwrapped'}))
    errors{end+1} = 'Phase output must be wrapped or unwrapped.';
end
selectionMethod = upper(strtrim(char(string(cfg.select_method))));
if ~any(strcmp(selectionMethod,{'PERCENT','DENSITY'}))
    errors{end+1} = 'Selection method must be PERCENT or DENSITY.';
elseif strcmp(selectionMethod,'PERCENT') && ...
        (~isfinite(cfg.percent_rand) || cfg.percent_rand < 0 || cfg.percent_rand > 100)
    errors{end+1} = 'Random percentage must be between 0 and 100.';
elseif strcmp(selectionMethod,'DENSITY') && ...
        (~isfinite(cfg.density_rand) || cfg.density_rand <= 0)
    errors{end+1} = 'Random density must be greater than zero.';
end
if any([cfg.weed_time_win cfg.unwrap_time_win cfg.scn_time_win] < 0)
    errors{end+1} = 'Temporal windows cannot be negative.';
end
if cfg.weed_standard_dev < 0 || cfg.weed_standard_dev > 100
    errors{end+1} = 'Weeding standard deviation must be between 0 and 100.';
elseif abs(cfg.weed_standard_dev * 10 - round(cfg.weed_standard_dev * 10)) > 1e-9
    errors{end+1} = 'Weeding standard deviation accepts at most one decimal place.';
end
if cfg.n_cores < 1 || mod(cfg.n_cores, 1) ~= 0
    errors{end+1} = 'Number of cores must be a positive integer.';
end
if isempty(strtrim(cfg.export_name))
    errors{end+1} = 'Export filename cannot be empty.';
end

if forStart
    runtime = phase_stamps_beta.inspectStaMPSInstallation( ...
        cfg.installation_folder,ispc);
    errors = [errors runtime.errors];
    if ~isfolder(cfg.project_path)
        errors{end+1} = ['PHASE project folder does not exist: ' cfg.project_path];
    end
else
    if ~isempty(cfg.installation_folder)
        runtime = phase_stamps_beta.inspectStaMPSInstallation( ...
            cfg.installation_folder,ispc);
        for k = 1:numel(runtime.errors)
            warnings{end+1} = runtime.errors{k}; %#ok<AGROW>
        end
    end
    if ~isempty(cfg.project_path) && ~isfolder(cfg.project_path)
        warnings{end+1} = 'The selected PHASE project folder does not currently exist.';
    end
end
end
