function [status, output] = runCommandHidden(app, command, label)
%RUNCOMMANDHIDDEN Run a shell command without opening a terminal window.
% Output is forwarded line by line to the PHASE StaMPS Run monitor.

if nargin < 3 || isempty(label), label = 'External command'; end
commandLine = buildCommandLine(command);
if ispc
    arguments = {'cmd.exe','/d','/q','/s','/c',['call ' commandLine]};
else
    arguments = {'/bin/bash','-lc',commandLine};
end

builder = java.lang.ProcessBuilder(javaStrings(arguments));
builder.redirectErrorStream(true);
environment = builder.environment();
% MATLAB setenv() updates the native process after the JVM has started.
% ProcessBuilder can otherwise retain the JVM's startup PATH and miss the
% StaMPS binary folders added immediately before this command.
if ispc
    pathKey = 'Path';
else
    pathKey = 'PATH';
end
environment.put(java.lang.String(pathKey), ...
    java.lang.String(getenv('PATH')));
environment.put(java.lang.String('PYTHONUNBUFFERED'),java.lang.String('1'));
environment.put(java.lang.String('PYTHONIOENCODING'),java.lang.String('utf-8'));

app.log(['Running ' char(string(label)) '…']);
try
    process = builder.start();
    reader = java.io.BufferedReader(java.io.InputStreamReader( ...
        process.getInputStream(),java.lang.String('UTF-8')));
catch ME
    error('PHASE_StaMPS:hiddenProcessStartFailed', ...
        'Could not start %s: %s',char(string(label)),ME.message);
end
readerCleanup = onCleanup(@() closeReader(reader)); %#ok<NASGU>

lines = {};
while process.isAlive() || reader.ready()
    while reader.ready()
        value = reader.readLine();
        if isempty(value), break; end
        line = char(value);
        lines{end+1,1} = line; %#ok<AGROW>
        if ~isempty(strtrim(line)), app.log(line); end
    end
    drawnow limitrate
    pause(0.04);
end
while reader.ready()
    value = reader.readLine();
    if isempty(value), break; end
    line = char(value);
    lines{end+1,1} = line; %#ok<AGROW>
    if ~isempty(strtrim(line)), app.log(line); end
end

status = double(process.waitFor());
output = strjoin(lines,newline);
if status == 0
    app.log([char(string(label)) ' completed.']);
end
end

function commandLine = buildCommandLine(command)
if iscell(command)
    values = cellfun(@(value) char(string(value)),command, ...
        'UniformOutput',false);
    if ispc
        values = cellfun(@quoteWindows,values,'UniformOutput',false);
    else
        values = cellfun(@quotePosix,values,'UniformOutput',false);
    end
    commandLine = strjoin(values,' ');
else
    commandLine = char(string(command));
end
end

function value = quoteWindows(value)
value = ['"' strrep(value,'"','""') '"'];
end

function value = quotePosix(value)
value = ['"' strrep(value,'"','\"') '"'];
end

function array = javaStrings(values)
array = javaArray('java.lang.String',numel(values));
for k = 1:numel(values)
    array(k) = java.lang.String(values{k});
end
end

function closeReader(reader)
try, reader.close(); catch, end
end
