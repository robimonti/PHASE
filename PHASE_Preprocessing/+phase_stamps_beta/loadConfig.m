function [cfg, info] = loadConfig(workDir)
%LOADCONFIG Load input_StaMPS.mat and migrate fields introduced later.

defaults = phase_stamps_beta.defaultConfig();
matPath = fullfile(workDir, 'input_StaMPS.mat');
info = struct('exists', false, 'migratedFields', {{}}, 'path', matPath);

if exist(matPath, 'file') ~= 2
    cfg = defaults;
    return
end

raw = load(matPath);
cfg = raw;
defaultNames = fieldnames(defaults);
for k = 1:numel(defaultNames)
    name = defaultNames{k};
    if ~isfield(cfg, name)
        if any(strcmp(name, {'weed_time_win','unwrap_time_win','scn_time_win'})) && ...
                isfield(cfg, 'time_span')
            cfg.(name) = cfg.time_span;
        else
            cfg.(name) = defaults.(name);
        end
        info.migratedFields{end+1} = name;
    end
end

cfg.stamps_first_step = char(string(cfg.stamps_first_step));
cfg.stamps_last_step = char(string(cfg.stamps_last_step));
cfg.master_date = char(string(cfg.master_date));
cfg.utc_time = char(string(cfg.utc_time));
cfg.export_name = char(string(cfg.export_name));
cfg.installation_folder = char(string(cfg.installation_folder));
cfg.project_path = char(string(cfg.project_path));
cfg.ph_output = char(string(cfg.ph_output));
info.exists = true;
end
