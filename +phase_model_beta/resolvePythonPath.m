function executable = resolvePythonPath(configured)
%RESOLVEPYTHONPATH Reuse the Python selected by the PHASE Windows installer.
%
% A saved input_model.mat remains authoritative. Before one exists, the
% standalone Model beta reads %APPDATA%\PHASE\python.txt, which is written by
% the installer and shared with the preprocessing and StaMPS modules.

if nargin < 1
    configured = '';
end
executable = stripOuterQuotes(strtrim(char(string(configured))));
if ~isempty(executable)
    return
end

if ispc
    appData = getenv('APPDATA');
    if ~isempty(appData)
        installerFile = fullfile(appData, 'PHASE', 'python.txt');
        if isfile(installerFile)
            executable = stripOuterQuotes(strtrim(fileread(installerFile)));
            if ~isempty(executable)
                return
            end
        end
    end
end

% Keep the same portable fallback used by the original Model app. The
% concrete interpreter can still be selected from the visible setup field.
executable = 'python3';
if ispc
    executable = 'python';
end
end

function value = stripOuterQuotes(value)
if numel(value) >= 2 && value(1) == '"' && value(end) == '"'
    value = value(2:end-1);
end
end
