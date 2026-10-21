function messages = ensureWindowsTriangleRuntime(stampsRoot)
%ENSUREWINDOWSTRIANGLERUNTIME Repair older managed installs before StaMPS.
% The upstream Windows binaries archive does not bundle MinGW DLLs. Keep
% Triangle's dependencies beside triangle.exe, avoiding a global PATH edit.

messages = {};
if ~ispc
    return
end

triangleBin = fullfile(stampsRoot,'external','triangle','bin');
triangleExe = fullfile(triangleBin,'triangle.exe');
if ~isfile(triangleExe)
    return % inspectStaMPSInstallation reports the missing executable.
end

localAppData = getenv('LOCALAPPDATA');
if isempty(localAppData)
    error('PHASE_StaMPS:triangleRuntimeMissing', ...
        'LOCALAPPDATA is unavailable; cannot locate the PHASE Portable Git runtime.');
end
mingwBin = fullfile(localAppData,'PHASE','portable-git','mingw64','bin');
dlls = {'libstdc++-6.dll','libgcc_s_seh-1.dll','libwinpthread-1.dll'};
for k = 1:numel(dlls)
    destination = fullfile(triangleBin,dlls{k});
    if isfile(destination)
        continue
    end
    source = fullfile(mingwBin,dlls{k});
    if ~isfile(source)
        error('PHASE_StaMPS:triangleRuntimeMissing', ...
            ['Triangle requires %s. It was not found in the PHASE Portable Git ', ...
             'runtime (%s). Repair the PHASE installation before running StaMPS.'], ...
            dlls{k},mingwBin);
    end
    [ok,message] = copyfile(source,destination);
    if ~ok
        error('PHASE_StaMPS:triangleRuntimeCopyFailed', ...
            'Could not install Triangle dependency %s: %s',dlls{k},message);
    end
    messages{end+1} = ['Triangle runtime repaired: ' dlls{k}]; %#ok<AGROW>
end

% The Windows port may return a nonzero code for -h even when it works.
% Exercise the same -e operation used by StaMPS Step 8 on three points.
probeDir = tempname;
[ok,message] = mkdir(probeDir);
if ~ok
    error('PHASE_StaMPS:triangleProbeSetupFailed', ...
        'Could not create Triangle probe folder: %s',message);
end
cleanup = onCleanup(@() rmdir(probeDir,'s')); %#ok<NASGU>
nodeFile = fullfile(probeDir,'phase-probe.1.node');
fid = fopen(nodeFile,'w');
if fid < 0
    error('PHASE_StaMPS:triangleProbeSetupFailed', ...
        'Could not create Triangle probe input: %s',nodeFile);
end
fprintf(fid,'3 2 0 0\n1 0 0\n2 1 0\n3 0 1\n');
fclose(fid);
command = ['cd /d "' probeDir '" && "' triangleExe ...
    '" -e phase-probe.1.node > triangle-probe.log 2>&1'];
[status,~] = system(command);
edgeFile = fullfile(probeDir,'phase-probe.2.edge');
if status ~= 0 || ~isfile(edgeFile)
    error('PHASE_StaMPS:triangleRuntimeFailed', ...
        ['Triangle could not triangulate a three-point test after its ', ...
         'dependencies were installed (exit %d). Repair the PHASE installation.'],status);
end
end
