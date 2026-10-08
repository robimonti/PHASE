function app = PHASE_StaMPS_beta(workDir)
%PHASE_STAMPS_BETA Launch the text-based PHASE StaMPS interface.
%
%   PHASE_StaMPS_beta opens the interface for the current folder.
%   PHASE_StaMPS_beta(WORKDIR) uses an explicit ASC_*/DSC_* folder.
%
% The interface is HTML/CSS/JavaScript, while all filesystem and StaMPS
% operations remain in editable MATLAB .m files. No App Designer round-trip
% or document.xml patch is required.

if nargin < 1 || isempty(workDir)
    workDir = pwd;
end
workDir = char(string(workDir));
if ~isfolder(workDir)
    error('PHASE_StaMPS_beta:workDirMissing', ...
        'StaMPS processing folder does not exist: %s',workDir);
end
if isfile(fullfile(workDir,'phase-project.json'))
    [~,projectPaths] = phase_project.open(workDir);
    datasets = [dir(fullfile(projectPaths.stamps,'ASC_*')); ...
        dir(fullfile(projectPaths.stamps,'DSC_*')); ...
        dir(fullfile(projectPaths.stamps,'DES_*'))];
    datasets = datasets([datasets.isdir]);
    if isempty(datasets)
        error('PHASE_StaMPS_beta:datasetMissing', ...
            'No StaMPS dataset exists yet in %s.',projectPaths.stamps);
    elseif numel(datasets) == 1
        workDir = fullfile(projectPaths.stamps,datasets(1).name);
    else
        selected = uigetdir(projectPaths.stamps,'Select a StaMPS dataset');
        if isequal(selected,0), app = []; return; end
        selected = char(java.io.File(selected).getCanonicalPath());
        valid = cellfun(@(name) fullfile(projectPaths.stamps,name), ...
            {datasets.name},'UniformOutput',false);
        if ispc
            isValid = any(strcmpi(selected,valid));
        else
            isValid = any(strcmp(selected,valid));
        end
        if ~isValid
            error('PHASE_StaMPS_beta:invalidDatasetSelection', ...
                'Select one of the StaMPS datasets in %s.',projectPaths.stamps);
        end
        workDir = selected;
    end
end

invocationDir = fileparts(mfilename('fullpath'));
launcherDir = resolveLauncherDir(invocationDir,workDir);
addpath(launcherDir);

app = phase_stamps_beta.App(workDir, launcherDir);
if nargout == 0
    clear app
end
end

function launcherDir = resolveLauncherDir(invocationDir,workDir)
% A convenience copy of this launcher may live inside ASC_*/DSC_*. Resolve
% the canonical editable package/UI instead of expecting assets beside it.
candidates = { ...
    invocationDir, ...
    fullfile(invocationDir,'PHASE_Preprocessing'), ...
    fullfile(invocationDir,'..','PHASE_Preprocessing'), ...
    fullfile(workDir,'PHASE_Preprocessing'), ...
    fullfile(workDir,'..','PHASE_Preprocessing')};
packageFile = which('phase_stamps_beta.App');
if ~isempty(packageFile)
    candidates{end+1} = fileparts(fileparts(packageFile));
end

for k = 1:numel(candidates)
    candidate = char(java.io.File(candidates{k}).getCanonicalPath());
    if isfile(fullfile(candidate,'+phase_stamps_beta','App.m')) && ...
            isfile(fullfile(candidate,'phase_stamps_beta_ui','index.html'))
        launcherDir = candidate;
        return
    end
end

error('PHASE_StaMPS_beta:runtimeAssetsMissing', ...
    ['The editable StaMPS beta runtime could not be found. Keep the complete ', ...
     'PHASE_Preprocessing/+phase_stamps_beta and phase_stamps_beta_ui folders ', ...
     'in the PHASE installation; copying only PHASE_StaMPS_beta.m is optional.']);
end
