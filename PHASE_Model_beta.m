function app = PHASE_Model_beta()
%PHASE_MODEL_BETA Launch the editable standalone PHASE Model application.
%
% The modern controller uses the complete mechanically extracted backend in
% +phase_model_beta/LegacyEngine.m. Runtime use does not read a legacy app.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
previousDir = pwd;
restoreDir = onCleanup(@() restoreFolder(previousDir)); %#ok<NASGU>
cd(rootDir);
app = phase_model_beta.App(rootDir);
if nargout == 0, clear app; end
end

function restoreFolder(folder)
try
    if isfolder(folder), cd(folder); end
catch
end
end
