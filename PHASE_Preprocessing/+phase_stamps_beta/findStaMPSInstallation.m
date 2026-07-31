function folder = findStaMPSInstallation(cfg, workDir)
%FINDSTAMPSINSTALLATION Locate a complete StaMPS clone near PHASE/dataset.

candidates = {};
if isfield(cfg,'installation_folder') && ~isempty(cfg.installation_folder)
    candidates{end+1} = cfg.installation_folder;
end
if isfield(cfg,'project_path') && ~isempty(cfg.project_path)
    candidates{end+1} = fullfile(cfg.project_path,'StaMPS');
    candidates{end+1} = fullfile(cfg.project_path,'engine','StaMPS');
end

candidates = [candidates, { ...
    fullfile(workDir,'StaMPS'), ...
    fullfile(workDir,'..','StaMPS'), ...
    fullfile(workDir,'..','engine','StaMPS'), ...
    fullfile(workDir,'..','..','StaMPS'), ...
    fullfile(workDir,'..','..','engine','StaMPS')}];

folder = '';
visited = {};
for k = 1:numel(candidates)
    candidate = canonicalPath(candidates{k});
    if any(samePath(candidate,visited))
        continue
    end
    visited{end+1} = candidate; %#ok<AGROW>
    if isStaMPSRoot(candidate)
        folder = candidate;
        return
    end
end
end

function tf = isStaMPSRoot(folder)
tf = isfolder(folder) && ...
    isfile(fullfile(folder,'matlab','stamps.m')) && ...
    isfile(fullfile(folder,'matlab','setparm.m'));
end

function pathValue = canonicalPath(pathValue)
pathValue = char(string(pathValue));
try
    pathValue = char(java.io.File(pathValue).getCanonicalPath());
catch
end
end

function matches = samePath(candidate,visited)
matches = false(size(visited));
for k = 1:numel(visited)
    if ispc
        matches(k) = strcmpi(candidate,visited{k});
    else
        matches(k) = strcmp(candidate,visited{k});
    end
end
end
