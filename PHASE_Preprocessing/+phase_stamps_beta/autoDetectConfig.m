function [cfg, detected, messages] = autoDetectConfig(cfg, workDir)
%AUTODETECTCONFIG Read acquisition metadata from SNAP .diff.par files.

detected = {};
messages = {};

candidates = {};
if isfield(cfg, 'project_path') && ~isempty(cfg.project_path)
    candidates{end+1} = fullfile(cfg.project_path, 'PHASE_Preprocessing');
end
candidates{end+1} = fullfile(workDir, '..', 'PHASE_Preprocessing');
candidates{end+1} = fullfile(workDir, 'PHASE_Preprocessing');

preprocFolder = '';
for k = 1:numel(candidates)
    candidate = char(java.io.File(candidates{k}).getCanonicalPath());
    if isfolder(candidate)
        matches = dir(fullfile(candidate, 'INSAR_*'));
        matches = matches([matches.isdir]);
        if ~isempty(matches)
            preprocFolder = candidate;
            break
        end
    end
end

if isempty(preprocFolder)
    messages{end+1} = 'Automatic detection: no PHASE_Preprocessing/INSAR_* folder found.';
    return
end

insarDirs = dir(fullfile(preprocFolder, 'INSAR_*'));
insarDirs = insarDirs([insarDirs.isdir]);
names = {insarDirs.name};
expected = ['INSAR_' char(string(cfg.master_date))];
index = find(strcmp(names, expected), 1);
if isempty(index)
    [~, index] = max([insarDirs.datenum]);
end
insarDir = insarDirs(index);
folderName = insarDir.name;
folderDate = extractAfter(folderName, 'INSAR_');
if ~isempty(regexp(folderDate, '^\d{8}$', 'once'))
    cfg.master_date = char(folderDate);
    detected{end+1} = 'master_date';
end

diff0Path = fullfile(preprocFolder, folderName, 'diff0');
parFiles = dir(fullfile(diff0Path, '*.diff.par'));
if isempty(parFiles)
    messages{end+1} = ['Automatic detection: no .diff.par files found in ' diff0Path];
    return
end

allDates = datetime.empty(1, 0);
for k = 1:numel(parFiles)
    tokens = regexp(parFiles(k).name, '(\d{8})_(\d{8})', 'tokens');
    if ~isempty(tokens)
        allDates(end+1) = datetime(tokens{1}{1}, 'InputFormat', 'yyyyMMdd'); %#ok<AGROW>
        allDates(end+1) = datetime(tokens{1}{2}, 'InputFormat', 'yyyyMMdd'); %#ok<AGROW>
    end
end
if ~isempty(allDates)
    minDate = min(allDates);
    maxDate = max(allDates);
    cfg.year_0 = year(minDate);
    cfg.month_0 = month(minDate);
    cfg.day_0 = day(minDate);
    cfg.time_span = days(maxDate - minDate);
    detected = [detected {'initial_date','time_span'}];
end

parText = fileread(fullfile(diff0Path, parFiles(1).name));
headMatch = regexp(parText, 'heading:\s+([-+\d\.Ee]+)', 'tokens', 'once');
if ~isempty(headMatch)
    cfg.heading = str2double(headMatch{1});
    detected{end+1} = 'heading';
end
freqMatch = regexp(parText, 'radar_frequency:\s+([-+\d\.Ee]+)', 'tokens', 'once');
if ~isempty(freqMatch)
    frequency = str2double(freqMatch{1});
    if isfinite(frequency) && frequency > 0
        cfg.lambda = 299792458 / frequency;
        detected{end+1} = 'lambda';
    end
end
timeMatch = regexp(parText, 'center_time:\s+([-+\d\.Ee]+)', 'tokens', 'once');
if ~isempty(timeMatch)
    secondsOfDay = str2double(timeMatch{1});
    hourValue = floor(secondsOfDay / 3600);
    minuteValue = round(floor(rem(secondsOfDay, 3600) / 60) / 10) * 10;
    if minuteValue == 60
        hourValue = hourValue + 1;
        minuteValue = 0;
    end
    hourValue = mod(hourValue, 24);
    cfg.utc_time = sprintf('%02d:%02d', hourValue, minuteValue);
    detected{end+1} = 'utc_time';
end

messages{end+1} = sprintf('Automatic detection: parameters read from %s.', diff0Path);
end
