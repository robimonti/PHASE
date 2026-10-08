function [project, p] = create(projectRoot, name)
%CREATE Create a data-only PHASE project. Never overwrite an existing one.

if nargin < 2 || isempty(name)
    [~,name] = fileparts(char(string(projectRoot)));
end
name = strtrim(char(string(name)));
if isempty(name)
    error('PHASE:ProjectNameMissing','Project name cannot be empty.');
end
p = phase_project.paths(projectRoot);
if isfile(p.manifest)
    error('PHASE:ProjectExists','A PHASE project already exists at %s.',p.root);
end
if isfolder(p.root)
    existing = dir(p.root);
    existing = existing(~ismember({existing.name},{'.','..'}));
    if ~isempty(existing)
        error('PHASE:ProjectFolderNotEmpty', ...
            'Select a new or empty project folder: %s.',p.root);
    end
else
    [ok,message] = mkdir(p.root);
    if ~ok, error('PHASE:ProjectCreateFailed','%s',message); end
end

folders = {p.raw,p.aoi,p.preprocessing,p.stamps,p.exports,p.model, ...
    p.figures,p.gis,p.reports,p.logs,p.imports};
for k = 1:numel(folders)
    [ok,message] = mkdir(folders{k});
    if ~ok, error('PHASE:ProjectCreateFailed','%s',message); end
end

project = struct( ...
    'schemaVersion',1, ...
    'id',char(java.util.UUID.randomUUID()), ...
    'name',name, ...
    'createdAt',char(datetime('now','TimeZone','UTC', ...
        'Format','yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
    'layout','phase-project-v1');
phase_project.writeJson(p.manifest,project);
end
