function valueType = tsValueType(workDir, cfg)
%TSVALUETYPE Resolve the series actually exported by the last PSI run.
valueType = '';
if strcmp(char(string(cfg.ph_output)), 'wrapped')
    return
end

metadataPath = fullfile(workDir,'EXPORT', ...
    [char(string(cfg.export_name)) '_series.json']);
if isfile(metadataPath)
    try
        metadata = jsondecode(fileread(metadataPath));
        candidate = char(string(metadata.valueType));
        if any(strcmp(candidate,{'v-do','v-dao'})) && ...
                isfile(fullfile(workDir,['ps_plot_ts_' candidate '.mat']))
            valueType = candidate;
            return
        end
    catch
        % A legacy or incomplete metadata file falls back to the MAT files.
    end
end

correctionRequested = cfg.train_flag == 0 && ...
    strcmpi(char(string(cfg.subtr_tropo)),'y');
if correctionRequested
    candidates = {'v-dao','v-do'};
else
    candidates = {'v-do','v-dao'};
end
for k = 1:numel(candidates)
    if isfile(fullfile(workDir,['ps_plot_ts_' candidates{k} '.mat']))
        valueType = candidates{k};
        return
    end
end
end
