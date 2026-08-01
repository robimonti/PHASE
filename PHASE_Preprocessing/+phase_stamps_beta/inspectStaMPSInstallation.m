function report = inspectStaMPSInstallation(folder, requireWindowsBinaries)
%INSPECTSTAMPSINSTALLATION Verify the MATLAB and native StaMPS runtime.

if nargin < 2
    requireWindowsBinaries = ispc;
end
folder = char(string(folder));
report = struct('ok',false,'folder',folder,'errors',{{}}, ...
    'warnings',{{}},'missingMatlab',{{}},'missingBinaries',{{}});

if isempty(strtrim(folder)) || ~isfolder(folder)
    report.errors{end+1} = ['StaMPS installation folder does not exist: ' folder];
    return
end

requiredMatlab = { ...
    fullfile('matlab','stamps.m'), ...
    fullfile('matlab','setparm.m'), ...
    fullfile('matlab','ps_load_initial.m')};
for k = 1:numel(requiredMatlab)
    if ~isfile(fullfile(folder,requiredMatlab{k}))
        report.missingMatlab{end+1} = requiredMatlab{k}; %#ok<AGROW>
    end
end
if isempty(findNamedFile(folder,ternary(ispc,'mt_prep_snap.bat','mt_prep_snap')))
    report.missingMatlab{end+1} = ternary(ispc,'mt_prep_snap.bat','mt_prep_snap');
end
if ~isempty(report.missingMatlab)
    report.errors{end+1} = sprintf( ...
        'The selected folder is not a complete StaMPS runtime. Missing: %s.', ...
        strjoin(report.missingMatlab,', '));
end

if requireWindowsBinaries
    requiredBinaries = { ...
        fullfile('bin','calamp.exe'), ...
        fullfile('bin','cpxsum.exe'), ...
        fullfile('bin','pscphase.exe'), ...
        fullfile('bin','pscdem.exe'), ...
        fullfile('bin','psclonlat.exe'), ...
        fullfile('bin','selpsc_patch.exe'), ...
        fullfile('bin','selsbc_patch.exe'), ...
        fullfile('external','triangle','bin','triangle.exe'), ...
        fullfile('external','snaphu','bin','snaphu.exe')};
    for k = 1:numel(requiredBinaries)
        if ~isfile(fullfile(folder,requiredBinaries{k}))
            report.missingBinaries{end+1} = requiredBinaries{k}; %#ok<AGROW>
        end
    end
    if ~isempty(report.missingBinaries)
        report.errors{end+1} = sprintf( ...
            ['The StaMPS Windows runtime is missing %d mandatory executable(s): %s. ', ...
             'Run prepare-windows-runtime.ps1 from the PHASE folder.'], ...
            numel(report.missingBinaries),strjoin(report.missingBinaries,', '));
    end
end

report.ok = isempty(report.errors);
end

function pathValue = findNamedFile(root,name)
pathValue = '';
matches = dir(fullfile(root,'**',name));
matches = matches(~[matches.isdir]);
if ~isempty(matches)
    pathValue = fullfile(matches(1).folder,matches(1).name);
end
end

function out = ternary(condition,yesValue,noValue)
if condition, out = yesValue; else, out = noValue; end
end
