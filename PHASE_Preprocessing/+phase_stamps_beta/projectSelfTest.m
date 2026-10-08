function projectSelfTest()
%PROJECTSELFTEST Verify the PHASE 7 StaMPS path and export adapters.

temporaryRoot = tempname();
mkdir(temporaryRoot);
cleanup = onCleanup(@() rmdir(temporaryRoot,'s')); %#ok<NASGU>
projectRoot = fullfile(temporaryRoot,'project');
[~,p] = phase_project.create(projectRoot,'StaMPS test');
dataset = fullfile(p.stamps,'ASC_Jan20_Feb20');
mkdir(fullfile(dataset,'EXPORT'));
exportSource = fullfile(dataset,'EXPORT','points.csv');
fid = fopen(exportSource,'w');
assert(fid > 0);
fprintf(fid,'id,displacement\n1,0\n');
fclose(fid);
insar = fullfile(p.preprocessing,'INSAR_20200101');
mkdir(insar);

assert(strcmp(phase_project.findRoot(dataset),p.root));
assert(strcmp(phase_stamps_beta.exportPath(p.root,'20200101'),insar));
published = phase_stamps_beta.publishExports(dataset);
assert(strcmp(published,fullfile(p.exports,'ASC_Jan20_Feb20')));
assert(isfile(fullfile(published,'points.csv')));
assert(isfile(exportSource));

cfg = phase_stamps_beta.defaultConfig();
[detected,fields] = phase_stamps_beta.autoDetectConfig(cfg,dataset);
assert(strcmp(detected.project_path,p.root));
assert(ismember('project_path',fields));
fprintf('PHASE StaMPS project adapter self-test passed.\n');
end
