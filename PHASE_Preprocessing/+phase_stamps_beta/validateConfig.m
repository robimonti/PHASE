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
if any([cfg.weed_time_win cfg.unwrap_time_win cfg.scn_time_win] < 0)
    errors{end+1} = 'Temporal windows cannot be negative.';
end
if cfg.n_cores < 1 || mod(cfg.n_cores, 1) ~= 0
    errors{end+1} = 'Number of cores must be a positive integer.';
end
if isempty(strtrim(cfg.export_name))
    errors{end+1} = 'Export filename cannot be empty.';
end

if forStart
    if ~isfolder(cfg.installation_folder)
        errors{end+1} = ['StaMPS installation folder does not exist: ' cfg.installation_folder];
    end
    if ~isfolder(cfg.project_path)
        errors{end+1} = ['PHASE project folder does not exist: ' cfg.project_path];
    end
else
    if ~isempty(cfg.installation_folder) && ~isfolder(cfg.installation_folder)
        warnings{end+1} = 'The selected StaMPS installation folder does not currently exist.';
    end
    if ~isempty(cfg.project_path) && ~isfolder(cfg.project_path)
        warnings{end+1} = 'The selected PHASE project folder does not currently exist.';
    end
end
end
