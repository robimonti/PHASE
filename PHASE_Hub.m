function app = PHASE_Hub(projectRoot)
%PHASE_HUB Open the unified PHASE workspace.
%   PHASE_Hub() opens the project chooser inside the hub.
%   PHASE_Hub(PROJECTROOT) opens an existing PHASE 7 project.

installation = phase_project.installationRoot();
addpath(installation);
addpath(fullfile(installation,'PHASE_Preprocessing'));
if nargin < 1, projectRoot = ''; end
app = phase_hub.App(installation,projectRoot);
if nargout == 0, clear app; end
end
