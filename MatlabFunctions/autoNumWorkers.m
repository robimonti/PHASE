function N = autoNumWorkers(maxSeriesLength)

% autoNumWorkers Restitutes a reasonable number of worker based on CPU and RAM available
%
% Input:
%   maxSeriesLength       - scalar with the length of the longest time series
%
% Output:
%   N                     - number of maximum workers

% Number of logic cores
try
    nCores = feature('numcores');
catch
    nCores = feature('numCores');  % fallback (older versions)
end

% RAM available (GB)
freeGB = 4; % fallback
try
    if ispc
        mem = memory;
        freeGB = mem.MaxPossibleArrayBytes / 1e9;
    elseif ismac || isunix
        % system commands to estimate RAM
        [~, res] = system('vm_stat');
        pageSize = 4096;   % byte per page
        pagesFree = regexp(res, 'Pages free:\s+(\d+)', 'tokens', 'once');
        pagesInactive = regexp(res, 'Pages inactive:\s+(\d+)', 'tokens', 'once');
        if ~isempty(pagesFree)
            pagesFree = str2double(pagesFree{1});
        else
            pagesFree = 0;
        end
        if ~isempty(pagesInactive)
            pagesInactive = str2double(pagesInactive{1});
        else
            pagesInactive = 0;
        end
        freeBytes = (pagesFree + pagesInactive) * pageSize;
        freeGB = freeBytes / 1e9;
    end
catch
    % safe fallback 
    warning('RAM detection failed. Using fallback value.');
end

% RAM necessary per worker (safe estimate)
RAM_per_worker = 0.2 + 0.00005 * maxSeriesLength^2;

% Maximum number of sustainable workers
maxRAM_workers = floor(freeGB / RAM_per_worker);

% Limit to the available number of cores
N = min([nCores, maxRAM_workers, 8]); % 8 = limit hard-coded

% Ensure at least 1
N = max(N, 1);
end