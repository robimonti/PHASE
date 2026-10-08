function runtime = runtime(projectRoot)
%RUNTIME Single source of paths for installed code and project data.

[project,p] = phase_project.open(projectRoot);
installation = phase_project.installationRoot();
runtime = struct();
runtime.project = project;
runtime.paths = p;
runtime.installation = installation;
runtime.preprocessingCode = fullfile(installation,'PHASE_Preprocessing');
runtime.modelCode = fullfile(installation,'+phase_model_beta');
runtime.pythonScripts = fullfile(installation,'pythonScripts');
runtime.snapGraphs = fullfile(installation,'PHASE_Preprocessing', ...
    'snap2stamps','graphs');
runtime.geoSplinter = fullfile(installation,'geoSplinter');
end
