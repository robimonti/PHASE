function [status, output] = runCommandLive(engine, command, phase, ...
        rangeStart, rangeEnd, firstStep, lastStep)
%RUNCOMMANDLIVE Run a generated preprocessing script without a terminal.
% stdout and stderr are merged and forwarded line-by-line to the beta Run
% monitor. The generated BAT/SH logic remains the same as in the stable app.

arguments
    engine
    command (1,:) char
    phase (1,:) char
    rangeStart (1,1) double = 0
    rangeEnd (1,1) double = 100
    firstStep (1,1) double = 1
    lastStep (1,1) double = 6
end

command = stripTrailingBackgroundOperator(command);
if ispc
    processArguments = {'cmd.exe','/d','/q','/s','/c', ...
        ['call "' command '"']};
else
    processArguments = {'/bin/bash',command};
end

builder = java.lang.ProcessBuilder(javaStrings(processArguments));
builder.redirectErrorStream(true);
environment = builder.environment();
environment.put(java.lang.String('PYTHONUNBUFFERED'),java.lang.String('1'));
environment.put(java.lang.String('PYTHONIOENCODING'),java.lang.String('utf-8'));

notifyProgress(engine,rangeStart,phase,false);
engine.updateOutput(['Running ' phase '…']);
started = tic;
lastNotice = tic;
outputLines = {};

try
    process = builder.start();
    engine.ActiveProcess = process;
    reader = java.io.BufferedReader(java.io.InputStreamReader( ...
        process.getInputStream(),java.lang.String('UTF-8')));
catch ME
    error('PHASE:LiveProcessStartFailed', ...
        'Could not start %s: %s',phase,ME.message);
end
processCleanup = onCleanup(@() clearActiveProcess(engine)); %#ok<NASGU>
readerCleanup = onCleanup(@() closeReader(reader)); %#ok<NASGU>

currentProgress = rangeStart;
while process.isAlive() || reader.ready()
    [lines,reader] = drainReadyLines(reader);
    for k = 1:numel(lines)
        line = lines{k};
        outputLines{end+1,1} = line; %#ok<AGROW>
        if ~isempty(strtrim(line))
            engine.updateOutput(line);
        end
        step = parseStep(line,firstStep,lastStep);
        if ~isempty(step)
            stepCount = max(1,lastStep-firstStep+1);
            fraction = max(0,min(1,(step-firstStep)/stepCount));
            currentProgress = max(currentProgress, ...
                rangeStart + fraction*(rangeEnd-rangeStart));
            notifyProgress(engine,currentProgress, ...
                sprintf('%s · step %d of %d',phase,step,lastStep),false);
        end
    end

    if toc(lastNotice) >= 0.5
        notifyProgress(engine,currentProgress,phase,false);
        lastNotice = tic;
    end
    drawnow limitrate
    pause(0.05);
end

% Drain the final bytes after process termination.
[lines,reader] = drainReadyLines(reader);
for k = 1:numel(lines)
    line = lines{k};
    outputLines{end+1,1} = line; %#ok<AGROW>
    if ~isempty(strtrim(line)), engine.updateOutput(line); end
end

status = double(process.waitFor());
output = strjoin(outputLines,newline);
elapsed = toc(started);
if engine.StopFlag
    notifyProgress(engine,currentProgress,[phase ' force-stopped'],false);
    error('PHASE:ProcessingStopped', ...
        '%s was force-stopped after %s.',phase,durationText(elapsed));
end
if status ~= 0
    notifyProgress(engine,currentProgress,[phase ' failed'],false);
    error('PHASE:ExternalCommandFailed', ...
        '%s failed with exit code %d after %s.%s%s', ...
        phase,status,durationText(elapsed),newline,lastUsefulOutput(outputLines));
end

notifyProgress(engine,rangeEnd,[phase ' completed'],false);
engine.updateOutput(sprintf('%s completed in %s.',phase,durationText(elapsed)));
end

function array = javaStrings(values)
array = javaArray('java.lang.String',numel(values));
for k = 1:numel(values)
    array(k) = java.lang.String(values{k});
end
end

function command = stripTrailingBackgroundOperator(command)
command = strtrim(command);
if endsWith(command,'&')
    command = strtrim(command(1:end-1));
end
end

function [lines,reader] = drainReadyLines(reader)
lines = {};
while reader.ready()
    value = reader.readLine();
    if isempty(value), break; end
    lines{end+1,1} = char(value); %#ok<AGROW>
end
end

function step = parseStep(line,firstStep,lastStep)
step = [];
token = regexp(line, ...
    '#{3,}\s*STEP\s+([1-6])(?:\s*\([^)]*\))?\s*#{3,}', ...
    'tokens','once','ignorecase');
if isempty(token), return; end
candidate = str2double(token{1});
if isfinite(candidate) && candidate >= firstStep && candidate <= lastStep
    step = candidate;
end
end

function notifyProgress(engine,percentage,phase,indeterminate)
if isempty(engine.ExternalProgressCallback), return; end
try
    engine.ExternalProgressCallback(struct( ...
        'percentage',max(0,min(100,double(percentage))), ...
        'phase',char(string(phase)), ...
        'indeterminate',logical(indeterminate)));
catch callbackError
    warning('PHASE:PreprocessingBetaProgressCallback', ...
        'Could not forward processing progress: %s',callbackError.message);
end
end

function text = durationText(seconds)
seconds = max(0,round(seconds));
hours = floor(seconds/3600);
minutes = floor(mod(seconds,3600)/60);
remaining = mod(seconds,60);
if hours > 0
    text = sprintf('%dh %02dm %02ds',hours,minutes,remaining);
elseif minutes > 0
    text = sprintf('%dm %02ds',minutes,remaining);
else
    text = sprintf('%ds',remaining);
end
end

function text = lastUsefulOutput(lines)
useful = lines(~cellfun(@(line) isempty(strtrim(line)),lines));
if isempty(useful)
    text = 'No diagnostic output was produced.';
else
    useful = useful(max(1,end-7):end);
    text = strjoin(useful,newline);
end
end

function closeReader(reader)
try, reader.close(); catch, end
end

function clearActiveProcess(engine)
try, engine.ActiveProcess = []; catch, end
end
