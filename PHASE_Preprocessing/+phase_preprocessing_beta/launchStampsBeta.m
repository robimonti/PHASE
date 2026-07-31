function app = launchStampsBeta(launcherFile, workDir)
%LAUNCHSTAMPSBETA Invoke the canonical launcher even if a stale copy exists
% inside the ASC_*/DSC_* dataset folder.

launcherFile = char(string(launcherFile));
workDir = char(string(workDir));
if ~isfile(launcherFile)
    error('PHASE_Preprocessing_beta:stampsLauncherMissing', ...
        'PHASE_StaMPS_beta launcher does not exist: %s',launcherFile);
end

launcherDir = fileparts(char(java.io.File(launcherFile).getCanonicalPath()));
previousDir = pwd;
restoreDir = onCleanup(@() restoreFolder(previousDir)); %#ok<NASGU>
addpath(launcherDir,'-begin');
cd(launcherDir);
app = PHASE_StaMPS_beta(workDir);
end

function restoreFolder(folder)
try
    if isfolder(folder), cd(folder); end
catch
end
end
