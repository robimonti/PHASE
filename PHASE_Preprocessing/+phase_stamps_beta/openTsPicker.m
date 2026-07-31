function pickerContainer = openTsPicker(workDir, cfg, parentContainer)
%OPENTSPICKER Load the native TS picker inside the PHASE StaMPS window.

if nargin < 3 || isempty(parentContainer) || ~isvalid(parentContainer)
    error('PHASE_StaMPS_beta:tsContainerMissing', ...
        'The integrated TS Points container is not available.');
end
pickerContainer = parentContainer;
if ~isfolder(workDir)
    error('PHASE_StaMPS_beta:workDirMissing', ...
        'Cannot resolve the StaMPS processing folder: %s', workDir);
end

if isfolder(cfg.installation_folder)
    addpath(fullfile(cfg.installation_folder, 'matlab'));
    addpath(fullfile(cfg.installation_folder, 'matlab_compat'));
end

valueType = 'v-do';
matPath = fullfile(workDir, ['ps_plot_ts_' valueType '.mat']);
if exist(matPath, 'file') ~= 2
    previous = pwd;
    cleanup = onCleanup(@() cd(previous));
    figuresBefore = findall(0, 'Type', 'figure');
    cd(workDir);
    evalc("ps_plot('" + valueType + "','ts',1)");
    figuresAfter = findall(0, 'Type', 'figure');
    newFigures = setdiff(figuresAfter, figuresBefore);
    for k = 1:numel(newFigures)
        try
            if isvalid(newFigures(k))
                delete(newFigures(k));
            end
        catch
        end
    end
    clear cleanup
    if exist(matPath, 'file') ~= 2
        error('PHASE_StaMPS_beta:tsDataMissing', ...
            'ps_plot finished but did not produce %s. Run StaMPS Step 7 first.', matPath);
    end
end

try
    delete(parentContainer.Children);
    ts_export_picker(workDir, parentContainer, valueType, ...
        char(string(cfg.export_name)));
catch ME
    rethrow(ME)
end
end
