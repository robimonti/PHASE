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
    p.figures,p.gis,p.reports,p.logs};
for k = 1:numel(folders)
    [ok,message] = mkdir(folders{k});
    if ~ok, error('PHASE:ProjectCreateFailed','%s',message); end
end

project = struct( ...
    'schemaVersion',3, ...
    'id',char(java.util.UUID.randomUUID()), ...
    'name',name, ...
    'createdAt',char(datetime('now','TimeZone','UTC', ...
        'Format','yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
    'layout','phase-project-v3');
phase_project.writeJson(p.manifest,project);
guide = fullfile(p.root,'README_PROJECT.txt');
fid = fopen(guide,'w');
if fid < 0
    error('PHASE:ProjectGuideFailed','Cannot create %s.',guide);
end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,['PHASE PROJECT FOLDERS\n\n', ...
    '01_INPUT               Your optional AOI files and SAR archive.\n', ...
    '02_PROCESSING_INTERNAL  PHASE working data; do not rename or edit during a run.\n', ...
    '03_RESULTS             Final PSI time series, models, figures, GIS and reports.\n', ...
    '04_LOGS                Diagnostic logs and optional legacy-import catalog.\n\n', ...
    'Use the PHASE app to import scenes, run processing and open results.\n']);
end
