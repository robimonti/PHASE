function [project, report] = PHASE_ImportLegacy(legacyRoot, projectRoot, mode)
%PHASE_IMPORTLEGACY Assisted import of final results from an older workspace.
%   PHASE_ImportLegacy prompts for folders and copies final results.
%   PHASE_ImportLegacy(OLD,NEW,'reference') catalogs them without copying.

if nargin < 1 || isempty(legacyRoot)
    legacyRoot = uigetdir(pwd,'Select the old PHASE workspace');
    if isequal(legacyRoot,0), project = []; report = []; return; end
end
if nargin < 2 || isempty(projectRoot)
    projectRoot = uigetdir(pwd,'Select an empty folder for the new project');
    if isequal(projectRoot,0), project = []; report = []; return; end
end
if nargin < 3, mode = 'copy'; end
inventory = phase_project.scanLegacy(legacyRoot);
fprintf('Found %d final files in %d Model runs and %d StaMPS datasets (%.2f GB).\n', ...
    numel(inventory.items),inventory.modelRuns,inventory.stampsDatasets, ...
    inventory.totalBytes / 1e9);
[project,report] = phase_project.importLegacy(legacyRoot,projectRoot,'',mode);
fprintf('Project created: %s\n',char(string(project.name)));
projectPaths = phase_project.paths(projectRoot);
fprintf('Import report: %s\n',fullfile(projectPaths.imports,'legacy-import.json'));
end
