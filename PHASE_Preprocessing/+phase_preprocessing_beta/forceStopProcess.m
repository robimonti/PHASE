function stopped = forceStopProcess(process)
%FORCESTOPPROCESS Terminate an external command and its complete child tree.

stopped = false;
if isempty(process), return; end
try
    if ~process.isAlive(), return; end
catch
    return
end

if ispc
    % taskkill /T includes cmd.exe, Python and every SNAP gpt.exe child.
    pid = sprintf('%.0f',double(process.pid()));
    arguments = {'taskkill.exe','/PID',pid,'/T','/F'};
    builder = java.lang.ProcessBuilder(javaStrings(arguments));
    builder.redirectErrorStream(true);
    killer = builder.start();
    killer.waitFor();
else
    destroyDescendants(process);
    try, process.destroyForcibly(); catch, process.destroy(); end
end

started = tic;
while toc(started) < 5
    try
        if ~process.isAlive()
            stopped = true;
            return
        end
    catch
        stopped = true;
        return
    end
    pause(0.05);
end

try, process.destroyForcibly(); catch, end
try, stopped = ~process.isAlive(); catch, stopped = true; end
end

function destroyDescendants(process)
% Java 9+ ProcessHandle API; silently fall back to killing only the parent.
try
    iterator = process.toHandle().descendants().iterator();
    handles = {};
    while iterator.hasNext()
        handles{end+1} = iterator.next(); %#ok<AGROW>
    end
    for k = numel(handles):-1:1
        try, handles{k}.destroyForcibly(); catch, end
    end
catch
end
end

function array = javaStrings(values)
array = javaArray('java.lang.String',numel(values));
for k = 1:numel(values)
    array(k) = java.lang.String(values{k});
end
end
