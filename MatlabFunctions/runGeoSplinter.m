function [status, output] = runGeoSplinter(executable, jobFile, stopCheck)
%RUNGEOSPLINTER Run a geoSplinter job without cmd.exe or shell redirection.
%
% Paths are resolved against the PHASE engine root. Java ProcessBuilder
% connects the job file directly to standard input, avoiding the temporary
% batch files and current-folder assumptions used by the legacy application.

if nargin < 3 || isempty(stopCheck)
    stopCheck = @() [];
end

runtimeRoot = phase_model_beta.projectRoot();
executable = canonicalPath(executable,runtimeRoot);
if ispc && ~endsWith(lower(executable),'.exe')
    executable = [executable '.exe'];
end
jobFile = canonicalPath(jobFile,runtimeRoot);

if ~isfile(executable)
    error('PHASE_Model_beta:geoSplinterMissing', ...
        'geoSplinter executable was not found: %s',executable);
end
if ~isfile(jobFile)
    error('PHASE_Model_beta:geoSplinterJobMissing', ...
        'geoSplinter job file was not created: %s',jobFile);
end

builder = java.lang.ProcessBuilder(javaStrings({executable}));
builder.directory(java.io.File(runtimeRoot));
builder.redirectInput(java.io.File(jobFile));
builder.redirectErrorStream(true);
try
    process = builder.start();
catch ME
    error('PHASE_Model_beta:geoSplinterStartFailed', ...
        'Could not start geoSplinter: %s',ME.message);
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
        try
            stopCheck();
        catch ME
            try, process.destroyForcibly(); catch, end
            rethrow(ME);
        end
        continue
    end
    text = char(line);
    lines{end+1} = text; %#ok<AGROW>
    fprintf('%s\n',text);
    drawnow limitrate
end
process.waitFor();
reader.close();
status = double(process.exitValue());
output = strjoin(lines,newline);
end

function value = canonicalPath(value,runtimeRoot)
candidate = java.io.File(char(string(value)));
if ~candidate.isAbsolute()
    candidate = java.io.File(runtimeRoot,char(string(value)));
end
value = char(candidate.getCanonicalPath());
end

function array = javaStrings(values)
array = javaArray('java.lang.String',numel(values));
for k = 1:numel(values)
    array(k) = java.lang.String(values{k});
end
end
