function state = defaultDownloader(rootDir, polygon)
%DEFAULTDOWNLOADER JSON-safe state for the integrated Sentinel-1 downloader.

filters = struct();
filters.dataset = 'SENTINEL-1';
filters.processingLevel = {'SLC'};
filters.beamMode = {'IW'};
filters.polarization = {};
filters.flightDirection = {};
filters.subtype = {};
filters.startDate = '';
filters.endDate = '';
filters.pathStart = '';
filters.pathEnd = '';
filters.frameStart = '';
filters.frameEnd = '';
filters.groupID = '';
filters.samplingRate = '';
filters.samplingUnit = 'Month';

lastPath = fullfile(rootDir, 'downloadasf', 'last_download_request.json');
if isfile(lastPath)
    try
        last = jsondecode(fileread(lastPath));
        names = fieldnames(filters);
        for k = 1:numel(names)
            name = names{k};
            if isfield(last, name), filters.(name) = last.(name); end
        end
        if isfield(last, 'sampling')
            if isfield(last.sampling, 'rate'), filters.samplingRate = last.sampling.rate; end
            if isfield(last.sampling, 'unit') && ~isempty(last.sampling.unit)
                filters.samplingUnit = last.sampling.unit;
            end
        end
        if isfield(last, 'aoi')
            if isfield(last.aoi, 'coordinates'), polygon = coordinateMatrix(last.aoi.coordinates); end
            if isfield(last.aoi, 'corners'), polygon = coordinateMatrix(last.aoi.corners); end
        end
    catch
    end
end

[loggedIn, username] = loginState(rootDir);
state = struct('filters', filters, 'polygon', polygon, ...
    'results', {repmat(resultTemplate(), 0, 1)}, ...
    'count', 0, 'totalSizeGB', 0, 'recommended', '', ...
    'status', 'Draw or reuse an AOI, choose filters, then search ASF.', ...
    'busy', false, 'progress', 0, 'loggedIn', loggedIn, 'username', username);
end

function template = resultTemplate()
template = struct('id','','sceneName','','date','','time','', ...
    'pathNumber',0,'frameNumber',0,'direction','','sizeGB',0, ...
    'sizeBytes',0,'selected',false,'coordinates',zeros(0,2));
end

function [loggedIn, username] = loginState(rootDir)
loggedIn = false; username = '';
resultPath = fullfile(rootDir, 'downloadasf', 'login_result.json');
requestPath = fullfile(rootDir, 'downloadasf', 'login_request.json');
if ~isfile(resultPath) || ~isfile(requestPath), return; end
try
    result = jsondecode(fileread(resultPath));
    request = jsondecode(fileread(requestPath));
    loggedIn = isfield(result, 'status') && strcmp(string(result.status), "success") && ...
        isfield(request, 'username') && strlength(string(request.username)) > 0;
    if loggedIn, username = char(string(request.username)); end
catch
end
end

function value = coordinateMatrix(raw)
while iscell(raw) && ~isempty(raw), raw = raw{1}; end
value = squeeze(double(raw));
if size(value,2) ~= 2 && size(value,1) == 2, value = value.'; end
if size(value,2) ~= 2, value = zeros(0,2); end
end
