function messages = prepareRuntime(cfg)
%PREPARERUNTIME Configure StaMPS/TRAIN paths in the current MATLAB process.

report = phase_stamps_beta.inspectStaMPSInstallation( ...
    cfg.installation_folder,ispc);
if ~report.ok
    error('PHASE_StaMPS_beta:incompleteStaMPSRuntime','%s', ...
        strjoin(report.errors,newline));
end

stampsRoot = canonicalPath(cfg.installation_folder);
addpath(genpath(stampsRoot));

binPaths = { ...
    fullfile(stampsRoot,'bin'), ...
    fullfile(stampsRoot,'external','triangle','bin'), ...
    fullfile(stampsRoot,'external','snaphu','bin')};
currentPath = strsplit(getenv('PATH'),pathsep);
for k = 1:numel(binPaths)
    if isfolder(binPaths{k}) && ~any(pathMatches(binPaths{k},currentPath))
        setenv('PATH',[binPaths{k} pathsep getenv('PATH')]);
        currentPath{end+1} = binPaths{k}; %#ok<AGROW>
    end
end

entryPoints = {'stamps','setparm','ps_load_initial'};
unresolved = {};
for k = 1:numel(entryPoints)
    if isempty(which(entryPoints{k}))
        unresolved{end+1} = entryPoints{k}; %#ok<AGROW>
    end
end
if ~isempty(unresolved)
    error('PHASE_StaMPS_beta:unresolvedStaMPSRuntime', ...
        'StaMPS was found but MATLAB cannot resolve: %s.', ...
        strjoin(unresolved,', '));
end

messages = {['StaMPS runtime ready: ' stampsRoot]};
trainRoot = findTrainRoot(cfg,stampsRoot);
if ~isempty(trainRoot)
    trainMatlab = fullfile(trainRoot,'matlab');
    if isfolder(trainMatlab)
        addpath(genpath(trainMatlab));
        messages{end+1} = ['TRAIN runtime ready: ' trainRoot];
    end
elseif cfg.train_flag == 0 && strcmpi(strtrim(cfg.subtr_tropo),'y')
    error('PHASE_StaMPS_beta:trainRuntimeMissing', ...
        ['TRAIN correction is enabled, but no TRAIN clone was found beside ', ...
         'PHASE/StaMPS. Run prepare-windows-runtime.ps1 or disable TRAIN.']);
end
end

function root = findTrainRoot(cfg,stampsRoot)
candidates = {fullfile(fileparts(stampsRoot),'TRAIN')};
if isfield(cfg,'project_path') && ~isempty(cfg.project_path)
    candidates = [candidates, { ...
        fullfile(cfg.project_path,'TRAIN'), ...
        fullfile(cfg.project_path,'engine','TRAIN')}];
end
root = '';
for k = 1:numel(candidates)
    candidate = canonicalPath(candidates{k});
    if isfolder(fullfile(candidate,'matlab'))
        root = candidate;
        return
    end
end
end

function pathValue = canonicalPath(pathValue)
pathValue = char(string(pathValue));
try
    pathValue = char(java.io.File(pathValue).getCanonicalPath());
catch
end
end

function matches = pathMatches(candidate,paths)
matches = false(size(paths));
for k = 1:numel(paths)
    if ispc
        matches(k) = strcmpi(candidate,paths{k});
    else
        matches(k) = strcmp(candidate,paths{k});
    end
end
end
