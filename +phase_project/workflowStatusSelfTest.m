function workflowStatusSelfTest()
%WORKFLOWSTATUSSELFTEST Ensure progress follows products, not empty folders.

temporaryRoot = tempname();
mkdir(temporaryRoot);
cleanup = onCleanup(@() rmdir(temporaryRoot,'s')); %#ok<NASGU>
[~,p] = phase_project.create(fullfile(temporaryRoot,'project'),'Progress test');
status = phase_project.workflowStatus(p.root);
assert(~status.preprocessing && ~status.stamps && ~status.model);
assert(strcmp(status.recommended,'preprocessing'));

dataset = fullfile(p.stamps,'ASC_example');
mkdir(dataset);
insar = fullfile(p.preprocessing,'INSAR_20200101');
mkdir(fullfile(insar,'diff0'));
mkdir(fullfile(insar,'geo'));
status = phase_project.workflowStatus(p.root);
assert(~status.preprocessing); % Empty output directories are not enough.
writeTestFile(fullfile(insar,'diff0','pair.diff'));
writeTestFile(fullfile(insar,'geo','lat.ras'));
status = phase_project.workflowStatus(p.root);
assert(status.preprocessing && strcmp(status.recommended,'stamps'));

published = fullfile(p.exports,'ASC_example');
mkdir(published);
writeTestFile(fullfile(published,'points.csv'));
status = phase_project.workflowStatus(p.root);
assert(~status.stamps);
writeTestFile(fullfile(published,'points.xlsx'));
status = phase_project.workflowStatus(p.root);
assert(status.stamps && strcmp(status.recommended,'model'));

modelRun = fullfile(p.model,'output_001');
mkdir(fullfile(modelRun,'files','mat'));
writeTestFile(fullfile(modelRun,'report.xlsx'));
status = phase_project.workflowStatus(p.root);
assert(~status.model);
writeTestFile(fullfile(modelRun,'files','mat','PHASEresults.mat'));
status = phase_project.workflowStatus(p.root);
assert(status.model && isempty(status.recommended));
fprintf('PHASE workflow status self-test passed.\n');
end

function writeTestFile(pathValue)
fid = fopen(pathValue,'w');
assert(fid > 0);
fprintf(fid,'test');
fclose(fid);
end
