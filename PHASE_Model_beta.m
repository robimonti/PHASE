function app = PHASE_Model_beta(projectRoot)
%PHASE_MODEL_BETA Launch the editable standalone PHASE Model application.
%
% The modern controller uses the complete mechanically extracted backend in
% +phase_model_beta/LegacyEngine.m. Runtime use does not read a legacy app.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
previousDir = pwd;
restoreDir = onCleanup(@() restoreFolder(previousDir)); %#ok<NASGU>
cd(rootDir);
if nargin < 1 || isempty(projectRoot)
    app = phase_model_beta.App(rootDir);
else
    phase_project.open(projectRoot);
    app = phase_model_beta.App(rootDir,projectRoot);
end
if nargout == 0, clear app; end
end

function restoreFolder(folder)
try
    if isfolder(folder), cd(folder); end
catch
end
end
