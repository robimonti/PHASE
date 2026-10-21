function [cfg, detected, messages] = autoDetectConfig(cfg, workDir)
%AUTODETECTCONFIG Read acquisition metadata from SNAP .diff.par files.

detected = {};
messages = {};

installationFolder = phase_stamps_beta.findStaMPSInstallation(cfg,workDir);
if ~isempty(installationFolder) && ...
        (~isfield(cfg,'installation_folder') || ...
         ~samePath(cfg.installation_folder,installationFolder))
    cfg.installation_folder = installationFolder;
    detected{end+1} = 'installation_folder';
    messages{end+1} = ['Automatic detection: StaMPS runtime found at ' ...
        installationFolder '.'];
end

candidates = {};
modernRoot = phase_project.findRoot(workDir);
if ~isempty(modernRoot)
    modernPaths = phase_project.paths(modernRoot);
    candidates{end+1} = modernPaths.preprocessing;
    if ~isfield(cfg,'project_path') || ~samePath(cfg.project_path,modernRoot)
        cfg.project_path = modernRoot;
        detected{end+1} = 'project_path';
    end
end
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

if isempty(modernRoot)
    projectFolder = char(java.io.File(fileparts(preprocFolder)).getCanonicalPath());
else
    projectFolder = modernRoot;
end
if ~isfield(cfg, 'project_path') || ~samePath(cfg.project_path, projectFolder)
    cfg.project_path = projectFolder;
    detected{end+1} = 'project_path';
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
% The SNAP2StaMPS .diff.par center_time can disagree with the actual
% acquisition UTC (observed 08:37 there for a 17:14 Sentinel-1 master).
% TRAIN needs the satellite pass time, so prefer the master SNAP metadata.
masterTime = masterAcquisitionUtc(preprocFolder, cfg.master_date);
timeMatch = regexp(parText, 'center_time:\s+([-+\d\.Ee]+)', 'tokens', 'once');
if ~isempty(masterTime)
    cfg.utc_time = masterTime;
    detected{end+1} = 'utc_time';
    messages{end+1} = ['Automatic detection: satellite UTC read from the SNAP master product (' masterTime ').'];
    if ~isempty(timeMatch)
        parMinutes = str2double(timeMatch{1}) / 60;
        masterMinutes = str2double(masterTime(1:2)) * 60 + str2double(masterTime(4:5));
        if isfinite(parMinutes) && abs(parMinutes - masterMinutes) > 5
            messages{end+1} = 'Warning: .diff.par center_time disagrees with the master acquisition UTC; the master value was used.';
        end
    end
elseif ~isempty(timeMatch)
    secondsOfDay = str2double(timeMatch{1});
    if isfinite(secondsOfDay) && secondsOfDay >= 0 && secondsOfDay < 86400
        hourValue = floor(secondsOfDay / 3600);
        minuteValue = floor(rem(secondsOfDay, 3600) / 60);
        cfg.utc_time = sprintf('%02d:%02d', hourValue, minuteValue);
        detected{end+1} = 'utc_time';
        messages{end+1} = 'Warning: satellite UTC came from .diff.par because master SNAP metadata was unavailable; verify it before using TRAIN.';
    end
end

messages{end+1} = sprintf('Automatic detection: parameters read from %s.', diff0Path);
detected = unique(detected, 'stable');
end

function utc = masterAcquisitionUtc(preprocFolder, masterDate)
utc = '';
masterFolder = fullfile(preprocFolder,'master');
if ~isfolder(masterFolder)
    return
end
products = dir(fullfile(masterFolder,'*.dim'));
masterDate = char(string(masterDate));
for k = 1:numel(products)
    if ~contains(products(k).name,masterDate)
        continue
    end
    productText = fileread(fullfile(products(k).folder,products(k).name));
    token = regexp(productText, ...
        '<MDATTR name="first_line_time"[^>]*>\s*\d{2}-[A-Za-z]{3}-\d{4}\s+(\d{2}:\d{2}):\d{2}', ...
        'tokens','once');
    if ~isempty(token)
        utc = token{1};
        return
    end
end
end

function tf = samePath(left, right)
left = char(string(left));
right = char(string(right));
if isempty(left) || isempty(right)
    tf = false;
    return
end
try
    left = char(java.io.File(left).getCanonicalPath());
    right = char(java.io.File(right).getCanonicalPath());
catch
end
if ispc
    tf = strcmpi(left, right);
else
    tf = strcmp(left, right);
end
end
