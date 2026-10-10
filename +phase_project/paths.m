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
layout = 'phase-project-v3';
if isfile(p.manifest)
    try
        manifest = jsondecode(fileread(p.manifest));
        layout = char(string(manifest.layout));
    catch
        % open() will report the invalid manifest; keep path resolution safe.
    end
end
p.layout = layout;
if strcmp(layout,'phase-project-v1')
    % Never rename an existing project: its processing state stays usable.
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
elseif strcmp(layout,'phase-project-v2')
    p.input = fullfile(projectRoot,'00_INPUT');
    p.raw = fullfile(p.input,'SAR_ARCHIVE');
    p.aoi = fullfile(p.input,'AOI');
    p.processing = fullfile(projectRoot,'10_PROCESSING_INTERNAL');
    p.preprocessing = fullfile(p.processing,'PREPROCESSING');
    p.stamps = fullfile(p.processing,'STAMPS_PSI');
    p.results = fullfile(projectRoot,'20_RESULTS');
    p.exports = fullfile(p.results,'PSI_TIME_SERIES');
    p.model = fullfile(p.results,'DISPLACEMENT_MODELS');
    p.figures = fullfile(p.results,'FIGURES');
    p.gis = fullfile(p.results,'GIS');
    p.reports = fullfile(p.results,'REPORTS');
    p.logs = fullfile(projectRoot,'90_LOGS');
    p.imports = fullfile(projectRoot,'95_IMPORTED_RESULTS');
else
    p.input = fullfile(projectRoot,'01_INPUT');
    p.raw = fullfile(p.input,'SAR_ARCHIVE');
    p.aoi = fullfile(p.input,'AOI');
    p.processing = fullfile(projectRoot,'02_PROCESSING_INTERNAL');
    p.preprocessing = fullfile(p.processing,'PREPROCESSING');
    p.stamps = fullfile(p.processing,'STAMPS_PSI');
    p.results = fullfile(projectRoot,'03_RESULTS');
    p.exports = fullfile(p.results,'PSI_TIME_SERIES');
    p.model = fullfile(p.results,'DISPLACEMENT_MODELS');
    p.figures = fullfile(p.results,'FIGURES');
    p.gis = fullfile(p.results,'GIS');
    p.reports = fullfile(p.results,'REPORTS');
    p.logs = fullfile(projectRoot,'04_LOGS');
    % Legacy-import metadata is diagnostic state, not a fifth user folder.
    p.imports = p.logs;
end
end
