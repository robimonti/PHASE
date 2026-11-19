function pythonDependenciesSelfTest(python)
%PYTHONDEPENDENCIESSELFTEST Managed first launch, cached check and stale marker.
% Pass an empty virtual environment to exercise real package installation.
root=tempname;mkdir(root);cleanup=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
engine=fullfile(root,'engine');mkdir(engine);
source=fullfile(phase_project.installationRoot(),'PHASE_Preprocessing','snap2stamps','bin','phase_python_dependencies.py');
target=fullfile(engine,'PHASE_Preprocessing','snap2stamps','bin');mkdir(target);copyfile(source,target);
writeJson(fullfile(root,'install.json'),struct('version','v7.2.2'));
old=getenv('PHASE_PYTHON');envCleanup=onCleanup(@()setenv('PHASE_PYTHON',old)); %#ok<NASGU>
setenv('PHASE_PYTHON',python);
phase_hub.ensurePythonDependencies(engine,[]);
marker=fullfile(root,'python-dependencies.json');assert(isfile(marker));
data=jsondecode(fileread(marker));assert(strcmp(data.python,python));
before=fileread(marker);phase_hub.ensurePythonDependencies(engine,[]);assert(strcmp(before,fileread(marker)));
writeJson(marker,struct('version','older','python',python));
phase_hub.ensurePythonDependencies(engine,[]);
data=jsondecode(fileread(marker));assert(~strcmp(data.version,'older'));
fprintf('Managed dependency installation, cached check and release change passed.\n');
end

function writeJson(path,data)
fid=fopen(path,'w');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(data));
end
