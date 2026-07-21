function saveConfig(workDir, cfg)
%SAVECONFIG Persist a configuration using the legacy MAT variable layout.

matPath = fullfile(workDir, 'input_StaMPS.mat');
if ~isfolder(workDir)
    error('PHASE_StaMPS_beta:workDirMissing', ...
        'Processing folder does not exist: %s', workDir);
end

temporary = fullfile(workDir, 'input_StaMPS.tmp.mat');
if exist(temporary, 'file') == 2
    builtin('delete', temporary);
end
save(temporary, '-struct', 'cfg', '-mat');

% Avoid leaving a partially-written configuration if MATLAB is interrupted.
if exist(matPath, 'file') == 2
    backup = [matPath '.bak'];
    copyfile(matPath, backup, 'f');
end
movefile(temporary, matPath, 'f');
end
