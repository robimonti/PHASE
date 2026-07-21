function pickerFigure = openTsPicker(workDir, cfg, ownerFigure)
%OPENTSPICKER Open the existing StaMPS picker from the web-based beta.

if nargin < 3, ownerFigure = []; end
pickerFigure = [];
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
            if isvalid(newFigures(k)) && ...
                    (isempty(ownerFigure) || newFigures(k) ~= ownerFigure)
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

pickerFigure = uifigure('Name', 'PHASE · TS Points', ...
    'Position', [120 90 1220 780], 'Color', [0.035 0.055 0.10]);
container = uipanel(pickerFigure, 'Position', [0 0 1220 780], ...
    'BorderType', 'none');
try
    ts_export_picker(workDir, container, valueType, char(string(cfg.export_name)));
catch ME
    if isvalid(pickerFigure), delete(pickerFigure); end
    rethrow(ME)
end
end
