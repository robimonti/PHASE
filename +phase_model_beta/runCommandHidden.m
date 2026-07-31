function [status, output] = runCommandHidden(app, commandLine, label)
%RUNCOMMANDHIDDEN Execute report helpers without opening a terminal window.

if nargin < 3 || isempty(label)
    label = 'External command';
end
commandLine = char(string(commandLine));
if ispc
    arguments = {'cmd.exe','/d','/q','/s','/c',['call ' commandLine]};
else
    arguments = {'/bin/sh','-lc',commandLine};
end

builder = java.lang.ProcessBuilder(javaStrings(arguments));
builder.redirectErrorStream(true);
try
    process = builder.start();
catch ME
    error('PHASE_Model_beta:commandStartFailed', ...
        'Could not start %s: %s', char(string(label)), ME.message);
end

reader = java.io.BufferedReader( ...
    java.io.InputStreamReader(process.getInputStream()));
lines = {};
while true
    line = reader.readLine();
    if isempty(line)
        if ~process.isAlive(), break; end
        pause(0.02);
        drawnow limitrate
        if app.StopRequested
            try, process.destroyForcibly(); catch, end
            phase_model_beta.throwIfStopped(app);
        end
        continue
    end
    text = char(line);
    lines{end+1} = text; %#ok<AGROW>
    forwardLog(app,text);
    drawnow limitrate
    if app.StopRequested
        try, process.destroyForcibly(); catch, end
        phase_model_beta.throwIfStopped(app);
    end
end
process.waitFor();
status = double(process.exitValue());
output = strjoin(lines,newline);
if ~isempty(output)
    output = [output newline];
end
end

function array = javaStrings(values)
array = javaArray('java.lang.String',numel(values));
for k = 1:numel(values)
    array(k) = java.lang.String(values{k});
end
end

function forwardLog(app, message)
try
    callback = app.ExternalLogCallback;
    if ~isempty(callback)
        callback(message);
    else
        fprintf('%s\n',message);
    end
catch
    fprintf('%s\n',message);
end
end
