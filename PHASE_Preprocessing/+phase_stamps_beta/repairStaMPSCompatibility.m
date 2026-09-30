function report = repairStaMPSCompatibility(installationFolder)
%REPAIRSTAMPSCOMPATIBILITY Apply narrowly scoped fixes for known upstream bugs.
%
% StaMPS master currently contains a typo in ps_select.m (sfprintf rather
% than sprintf). It is reached only when gamma_stdev_reject > 0, which made
% normal configurations fail at Step 3 despite all PHASE inputs being valid.
% Keep a one-time backup and change only that exact token.

report = struct('changed', false, 'message', '');
psSelect = fullfile(char(string(installationFolder)), 'matlab', 'ps_select.m');
if ~isfile(psSelect)
    report.message = ['ps_select.m was not found at ' psSelect];
    return;
end

source = fileread(psSelect);
if ~contains(source, 'sfprintf(')
    return;
end

backup = [psSelect '.phase-original'];
if ~isfile(backup)
    [ok, message] = copyfile(psSelect, backup);
    if ~ok
        error('PHASE_StaMPS:compatibilityBackupFailed', ...
            'Could not back up %s before applying the StaMPS compatibility fix: %s', ...
            psSelect, message);
    end
end

fixed = strrep(source, 'sfprintf(', 'sprintf(');
fid = fopen(psSelect, 'w');
if fid == -1
    error('PHASE_StaMPS:compatibilityWriteFailed', ...
        'Could not write the StaMPS compatibility fix to %s.', psSelect);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', fixed);

clear ps_select
report.changed = true;
report.message = ['Applied StaMPS compatibility fix in ' psSelect ...
    ' (sfprintf -> sprintf; original saved as .phase-original).'];
end
