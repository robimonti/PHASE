function app = PHASE_Preprocessing_beta()
%PHASE_PREPROCESSING_BETA Launch the text-based PHASE preprocessing UI.
%
% Run this file from MATLAB. The beta has no runtime dependency on either
% legacy MLAPP; it uses a mechanically extracted, editable MATLAB engine.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
addpath(fullfile(rootDir, 'PHASE_Preprocessing'));
app = phase_preprocessing_beta.App(rootDir);
if nargout == 0, clear app; end
end
