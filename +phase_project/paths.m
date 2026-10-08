function p = paths(projectRoot)
%PATHS Resolve the PHASE project layout from one canonical root.

projectRoot = char(string(projectRoot));
if isempty(strtrim(projectRoot))
    error('PHASE:ProjectPathMissing','Provide a project folder.');
end
projectRoot = char(java.io.File(projectRoot).getCanonicalPath());
p = struct();
p.root = projectRoot;
p.manifest = fullfile(projectRoot,'phase-project.json');
p.input = fullfile(projectRoot,'input');
p.raw = fullfile(p.input,'raw');
p.aoi = fullfile(p.input,'aoi');
p.processing = fullfile(projectRoot,'processing');
p.preprocessing = fullfile(p.processing,'preprocessing');
p.stamps = fullfile(p.processing,'stamps');
p.results = fullfile(projectRoot,'results');
p.exports = fullfile(p.results,'exports');
p.model = fullfile(p.results,'model');
p.figures = fullfile(p.results,'figures');
p.gis = fullfile(p.results,'gis');
p.reports = fullfile(p.results,'reports');
p.logs = fullfile(projectRoot,'logs');
p.imports = fullfile(projectRoot,'imports');
end
