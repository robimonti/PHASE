function app = PHASE_Preprocessing_beta(projectRoot)
%PHASE_PREPROCESSING_BETA Launch the text-based PHASE preprocessing UI.
%
% With a PHASE 7 project root, code stays in the installation while data
% and generated SNAP/StaMPS files remain in the selected project.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
addpath(fullfile(rootDir, 'PHASE_Preprocessing'));
if nargin < 1 || isempty(projectRoot)
    projectRoot = rootDir;
else
    projectRoot = char(java.io.File(char(string(projectRoot))).getCanonicalPath());
    phase_project.open(projectRoot);
end
app = phase_preprocessing_beta.App(rootDir,projectRoot);
if nargout == 0, clear app; end
end
