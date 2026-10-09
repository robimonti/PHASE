function selfTest()
%SELFTEST Exercise project creation, validation, scanning and import.

testRoot = tempname();
mkdir(testRoot);
cleanup = onCleanup(@() rmdir(testRoot,'s')); %#ok<NASGU>
legacy = fullfile(testRoot,'legacy');
mkdir(fullfile(legacy,'output_001'));
mkdir(fullfile(legacy,'ASC_Jan20_Feb20','EXPORT'));
mkdir(fullfile(legacy,'PHASE_Preprocessing','slaves'));
writeFixture(fullfile(legacy,'output_001','report.xlsx'));
writeFixture(fullfile(legacy,'ASC_Jan20_Feb20','EXPORT','points.csv'));
writeFixture(fullfile(legacy,'PHASE_Preprocessing','slaves','raw.zip'));

inventory = phase_project.scanLegacy(legacy);
assert(numel(inventory.items) == 2);
assert(inventory.modelRuns == 1 && inventory.stampsDatasets == 1);
assert(~contains(strjoin({inventory.items.source},'|'),'raw.zip'));

destination = fullfile(testRoot,'new-project');
[project,report] = phase_project.importLegacy(legacy,destination,'Test','copy');
[loaded,p] = phase_project.open(destination);
assert(strcmp(project.id,loaded.id));
assert(strcmp(loaded.name,'Test'));
assert(strcmp(loaded.layout,'phase-project-v2'));
assert(isfile(fullfile(p.root,'README_PROJECT.txt')));
assert(contains(p.results,'20_RESULTS'));
assert(numel(report.items) == 2);
assert(isfile(fullfile(p.model,'output_001','report.xlsx')));
assert(isfile(fullfile(p.exports,'ASC_Jan20_Feb20','points.csv')));
assert(isfile(fullfile(p.imports,'legacy-import.json')));
assert(isfile(fullfile(legacy,'output_001','report.xlsx')));
assert(strcmp(phase_project.runtime(destination).paths.root,p.root));
copiedResults = phase_project.importedResults(destination);
assert(numel(copiedResults) == 2 && all([copiedResults.exists]));
assert(all([copiedResults.copied]));

referenceDestination = fullfile(testRoot,'reference-project');
[~,referenceReport] = phase_project.importLegacy( ...
    legacy,referenceDestination,'Reference','reference');
assert(strcmp(referenceReport.mode,'reference'));
referencePaths = phase_project.paths(referenceDestination);
assert(~isfile(fullfile(referencePaths.model,'output_001','report.xlsx')));
referencedResults = phase_project.importedResults(referenceDestination);
assert(numel(referencedResults) == 2 && all([referencedResults.exists]));
assert(~any([referencedResults.copied]));

legacyProjectRoot = fullfile(testRoot,'existing-v1-project');
mkdir(legacyProjectRoot);
phase_project.writeJson(fullfile(legacyProjectRoot,'phase-project.json'), ...
    struct('schemaVersion',1,'id','existing','name','Existing', ...
    'layout','phase-project-v1'));
[~,legacyPaths] = phase_project.open(legacyProjectRoot);
assert(strcmp(legacyPaths.processing,fullfile(legacyProjectRoot,'processing')));
assert(strcmp(legacyPaths.results,fullfile(legacyProjectRoot,'results')));

installation = phase_project.installationRoot();
uiCache = fullfile(testRoot,'ui-cache');
sources = { ...
    fullfile(installation,'phase_model_beta_ui'), ...
    fullfile(installation,'PHASE_Preprocessing','phase_preprocessing_beta_ui'), ...
    fullfile(installation,'PHASE_Preprocessing','phase_stamps_beta_ui')};
modules = {'model','preprocessing','stamps'};
for k = 1:numel(modules)
    uiDir = phase_project.uiRuntime(sources{k},modules{k},uiCache);
    assert(isfile(fullfile(uiDir,'index.html')));
    assert(isfolder(fullfile(uiDir,'runtime_logs')));
    assert(isfolder(fullfile(uiDir,'map_tiles')));
end
fprintf('PHASE project self-test passed.\n');
end

function writeFixture(pathValue)
fid = fopen(pathValue,'w');
assert(fid > 0);
fprintf(fid,'fixture\n');
fclose(fid);
end
