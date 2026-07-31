function [executable, info] = resolvePython(configured)
%RESOLVEPYTHON Resolve a real Python 3 executable for PHASE backends.
%
% On Windows the PHASE installer records the selected interpreter in
% %APPDATA%\PHASE\python.txt. Prefer that file over an ambiguous `python`
% command, which may still point to Python 2 on machines running StaMPS.

if nargin < 1, configured = ''; end
configured = stripOuterQuotes(strtrim(char(string(configured))));
candidates = {};
sources = {};

isGeneric = any(strcmpi(configured, ...
    {'','python','python2','python2.7','python3','python3.11'}));
if ~isGeneric
    [candidates, sources] = addCandidate(candidates, sources, configured, 'configured value');
end

if ispc
    appData = getenv('APPDATA');
    if ~isempty(appData)
        installerFile = fullfile(appData, 'PHASE', 'python.txt');
        if isfile(installerFile)
            installerPython = stripOuterQuotes(strtrim(fileread(installerFile)));
            [candidates, sources] = addCandidate(candidates, sources, ...
                installerPython, '%APPDATA%\PHASE\python.txt');
        end
    end
    [candidates, sources] = addCandidate(candidates, sources, 'py -3', 'Windows Python launcher');
end

if ~isempty(configured)
    [candidates, sources] = addCandidate(candidates, sources, configured, 'configured value');
end
[candidates, sources] = addCandidate(candidates, sources, 'python3', 'python3 on PATH');
[candidates, sources] = addCandidate(candidates, sources, 'python', 'python on PATH');

attempts = {};
probeCode = [ ...
    'import sys; print(''PHASE_PYTHON=%d.%d|%s'' % ' ...
    '(sys.version_info[0], sys.version_info[1], sys.executable))'];
for k = 1:numel(candidates)
    command = candidateCommand(candidates{k});
    [status, output] = system([command ' -c ' quoteArgument(probeCode)]);
    token = regexp(output, ...
        'PHASE_PYTHON=(\d+)\.(\d+)\|([^\r\n]+)', 'tokens', 'once');
    if status == 0 && ~isempty(token)
        major = str2double(token{1});
        minor = str2double(token{2});
        resolved = stripOuterQuotes(strtrim(token{3}));
        if major == 3 && minor >= 8
            % pythonw.exe suppresses/loses stdout in GUI applications. The
            % preprocessing beta captures a sibling python.exe invisibly, so
            % SNAP/Python messages can be streamed into Run monitor.
            if ispc && endsWith(lower(resolved),'pythonw.exe')
                consolePython = fullfile(fileparts(resolved),'python.exe');
                if isfile(consolePython), resolved = consolePython; end
            end
            executable = resolved;
            info = struct( ...
                'source', sources{k}, ...
                'version', sprintf('%d.%d', major, minor), ...
                'configured', configured, ...
                'changed', ~sameExecutable(configured, resolved), ...
                'attempts', {attempts});
            return
        end
        attempts{end+1} = sprintf('%s resolved to unsupported Python %d.%d', ...
            candidates{k}, major, minor); %#ok<AGROW>
    else
        attempts{end+1} = sprintf('%s was unavailable', candidates{k}); %#ok<AGROW>
    end
end

error('PHASE_Preprocessing_beta:python3NotFound', ...
    ['PHASE requires Python 3.8 or newer, but no compatible interpreter was found. ' ...
     'The configured value was "%s". Run the PHASE installer again or set ' ...
     'the Setup > Python field to the full path of python.exe. Tried: %s.'], ...
    configured, strjoin(attempts, '; '));
end

function [candidates, sources] = addCandidate(candidates, sources, value, source)
value = stripOuterQuotes(strtrim(char(string(value))));
if isempty(value), return; end
if any(strcmpi(value, candidates)), return; end
candidates{end+1} = value;
sources{end+1} = source;
end

function command = candidateCommand(candidate)
standardCommands = {'python','python2','python2.7','python3','python3.11','py -3'};
if any(strcmpi(candidate, standardCommands))
    command = candidate;
else
    command = quoteArgument(candidate);
end
end

function value = quoteArgument(raw)
value = char(string(raw));
value = strrep(value, '"', '""');
value = ['"' value '"'];
end

function value = stripOuterQuotes(value)
if numel(value) >= 2 && value(1) == '"' && value(end) == '"'
    value = value(2:end-1);
end
end

function tf = sameExecutable(left, right)
left = stripOuterQuotes(strtrim(char(string(left))));
right = stripOuterQuotes(strtrim(char(string(right))));
if ispc, tf = strcmpi(left, right); else, tf = strcmp(left, right); end
end
