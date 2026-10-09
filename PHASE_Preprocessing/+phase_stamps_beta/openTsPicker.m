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

valueType = phase_stamps_beta.tsValueType(workDir,cfg);
if isempty(valueType)
    if strcmp(char(string(cfg.ph_output)),'wrapped')
        error('PHASE_StaMPS_beta:wrappedNoPicker', ...
            'TS Points is for unwrapped displacement. Wrapped phase is already exported in EXPORT.');
    end
    error('PHASE_StaMPS_beta:tsDataMissing', ...
        'No displacement time series is available yet. Complete StaMPS Step 7 first.');
end
try
    delete(parentContainer.Children);
    ts_export_picker(workDir, parentContainer, valueType, ...
        char(string(cfg.export_name)));
catch ME
    rethrow(ME)
end
end
