function context = sentinelUpdateContext(rootDir)
%SENTINELUPDATECONTEXT Find the reference and dates used by stack update.

projectFolder = fullfile(rootDir,'PHASE_Preprocessing');
slavesFolder = fullfile(projectFolder,'slaves');
reference = ''; dates = {};

masterFiles = dir(fullfile(projectFolder,'master','*.zip'));
slaveFiles = dir(fullfile(slavesFolder,'**','*.zip'));
allFiles = [masterFiles; slaveFiles];
for k = 1:numel(allFiles)
    pathValue = fullfile(allFiles(k).folder,allFiles(k).name);
    token = regexp(allFiles(k).name,'_(\d{8})T\d{6}_','tokens','once');
    if ~isempty(token), dates{end+1} = token{1}; end %#ok<AGROW>
    if isempty(reference), reference = pathValue; end
end

directories = dir(slavesFolder);
for k = 1:numel(directories)
    if directories(k).isdir && ~isempty(regexp(directories(k).name,'^\d{8}$','once'))
        dates{end+1} = directories(k).name; %#ok<AGROW>
    end
end
dates = unique(dates);
if isempty(dates), latest = ''; else
    values = datetime(dates,'InputFormat','yyyyMMdd'); latest = datestr(max(values),'yyyymmdd');
end
available = ~isempty(latest) && ~isempty(reference);
if isempty(latest), message = 'No local Sentinel-1 acquisition date is available.';
elseif isempty(reference), message = 'A reference Sentinel-1 ZIP could not be found.';
else, message = sprintf('Search newer acquisitions after %s using %s.',latest,basename(reference)); end
context = struct('available',available,'latestDate',latest,'referenceZip',reference, ...
    'localDates',{dates},'message',message);
end

function value = basename(pathValue)
[name,base,ext] = fileparts(pathValue); %#ok<ASGLU>
value = [base ext];
end
