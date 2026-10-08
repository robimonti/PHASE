function command = scriptCommand(python, scriptName, configPath, installRoot)
%SCRIPTCOMMAND Invoke installed SNAP Python code with project configuration.

scriptName = char(string(scriptName));
if ~isempty(regexp(scriptName, '[\\/]', 'once')) || ...
        ~endsWith(scriptName, '.py')
    error('PHASE:InvalidSnapScript', 'Expected a SNAP Python script filename.');
end
scriptPath = fullfile(char(string(installRoot)), 'PHASE_Preprocessing', ...
    'snap2stamps', 'bin', scriptName);
if ~isfile(scriptPath)
    error('PHASE:SnapScriptMissing', 'SNAP Python script not found: %s', scriptPath);
end
values = {char(string(python)), scriptPath, char(string(configPath))};
if any(cellfun(@(value) isempty(value) || contains(value, newline), values))
    error('PHASE:InvalidSnapCommand', 'Python and script paths must be nonempty single lines.');
end
if ispc
    if any(cellfun(@(value) contains(value, '"'), values))
        error('PHASE:InvalidSnapCommand', 'Double quotes are not supported in Windows paths.');
    end
    quoted = cellfun(@(value) ['"' value '"'], values, 'UniformOutput', false);
else
    quoted = cellfun(@(value) ['''' strrep(value, '''', '''"''"''') ''''], ...
        values, 'UniformOutput', false);
end
command = strjoin(quoted, ' ');
if ispc
    command = [command ' || exit /b 1'];
end
end
