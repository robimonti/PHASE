function projectPathsSelfTest()
%PROJECTPATHSSELFTEST Verify data helpers in an isolated PHASE 7 project.

temporaryRoot = tempname();
mkdir(temporaryRoot);
cleanup = onCleanup(@() rmdir(temporaryRoot,'s')); %#ok<NASGU>
projectRoot = fullfile(temporaryRoot,'project');
[~,p] = phase_project.create(projectRoot,'Preprocessing test');
assert(strcmp(phase_preprocessing_beta.dataFolder(p.root),p.preprocessing));

cfg = phase_preprocessing_beta.defaultConfig();
cfg.master_date = '20200101';
phase_preprocessing_beta.saveConfig(p.root,cfg);
[loaded,info] = phase_preprocessing_beta.loadConfig(p.root);
assert(info.exists && strcmp(loaded.master_date,cfg.master_date));
assert(isfile(fullfile(p.preprocessing,'input_preprocessing.mat')));

slaves = fullfile(p.preprocessing,'slaves');
mkdir(slaves);
fileName = 'S1A_IW_SLC__1SDV_20200101T000000_20200101T000025_test.zip';
fid = fopen(fullfile(slaves,fileName),'w');
assert(fid > 0);
fclose(fid);
inventory = phase_preprocessing_beta.scanSlaves(p.root);
assert(numel(inventory) == 1 && strcmp(inventory(1).name,fileName));
context = phase_preprocessing_beta.sentinelUpdateContext(p.root);
assert(strcmp(context.latestDate,'20200101'));
assert(isempty(phase_preprocessing_beta.collectFootprints(p.root)));

% A wrapper in a project with spaces must run installed Python source with
% an absolute project configuration path, independent of the current folder.
installRoot = fullfile(temporaryRoot,'installation with spaces');
scriptFolder = fullfile(installRoot,'PHASE_Preprocessing','snap2stamps','bin');
mkdir(scriptFolder);
scriptPath = fullfile(scriptFolder,'SEN_master_selector.py');
fid = fopen(scriptPath,'w');
assert(fid > 0);
fprintf(fid,'import pathlib, sys\nprint(pathlib.Path(sys.argv[1]).read_text().strip())\n');
fclose(fid);
configPath = fullfile(p.preprocessing,'snap2stamps','bin','project_master.conf');
mkdir(fileparts(configPath));
fid = fopen(configPath,'w');
assert(fid > 0);
fprintf(fid,'PROJECT_PATHS_OK\n');
fclose(fid);
python = phase_preprocessing_beta.resolvePython('');
command = phase_preprocessing_beta.scriptCommand(python, ...
    'SEN_master_selector.py',configPath,installRoot);
[status,output] = system(command);
assert(status == 0 && contains(output,'PROJECT_PATHS_OK'));
assert(strcmp(phase_preprocessing_beta.stampsFolder(p.root),p.stamps));
assert(isfolder(p.stamps));
fprintf('PHASE Preprocessing project path self-test passed.\n');
end
