function [project, paths] = PHASE_NewProject(projectRoot, name)
%PHASE_NEWPROJECT Create an empty data-only PHASE project.
%   PHASE_NewProject prompts for an empty destination folder.

if nargin < 1 || isempty(projectRoot)
    projectRoot = uigetdir(pwd,'Select an empty folder for the new PHASE project');
    if isequal(projectRoot,0), project = []; paths = []; return; end
end
if nargin < 2, name = ''; end
[project,paths] = phase_project.create(projectRoot,name);
fprintf('PHASE project created: %s\n',paths.root);
end
