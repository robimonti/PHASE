classdef App < handle
    %APP Modern controller around the complete text-extracted preprocessing engine.

    properties (SetAccess = private)
        UIFigure
        HTML
        Engine
        RootDir
        Config
        SavedConfig
        Logs = {}
        Status = 'idle'
        StatusDetail = 'Initialising preprocessing engine'
        IsDirty = false
        IsRunning = false
        StopRequested = false
        RunStartedAt = []
        RunProgress = struct('percentage',0,'phase','Ready', ...
            'elapsedSeconds',0,'etaSeconds',-1,'startedAt','', ...
            'indeterminate',false)
        MapBase
        MapFootprints
        MapPolygon
        SlaveFiles
        Downloader
        UpdateContext
        UpdateResults = struct([])
        UpdateStatus = 'Load the local stack context to search for newer acquisitions.'
        UpdateBusy = false
        MapTileBusy = false
        Transfer
        TransferProcess = []
        TransferTimer = []
        TransferDirectory = ''
        TransferLogPath = ''
    end

    methods
        function obj = App(rootDir)
            obj.RootDir = char(java.io.File(rootDir).getCanonicalPath());
            [loadedConfig, info] = phase_preprocessing_beta.loadConfig(obj.RootDir);
            obj.Config = loadedConfig;
            obj.SavedConfig = loadedConfig;
            pythonNotice = '';
            pythonChanged = false;
            try
                [resolvedPython, pythonInfo] = ...
                    phase_preprocessing_beta.resolvePython(obj.Config.python);
                obj.Config.python = resolvedPython;
                pythonChanged = pythonInfo.changed;
                pythonNotice = sprintf('Python %s selected from %s: %s', ...
                    pythonInfo.version, pythonInfo.source, resolvedPython);
            catch ME
                pythonNotice = ['Python 3 is not configured: ' ME.message];
            end
            gptNotice = '';
            gptChanged = false;
            try
                [resolvedGpt,gptInfo] = ...
                    phase_preprocessing_beta.resolveGpt(obj.Config.gptbin_path);
                obj.Config.gptbin_path = resolvedGpt;
                gptChanged = gptInfo.changed;
                if gptInfo.exists
                    gptNotice = ['SNAP gpt selected: ' resolvedGpt];
                else
                    gptNotice = ['SNAP gpt is not present at the configured path: ' resolvedGpt];
                end
            catch ME
                gptNotice = ['SNAP gpt is not configured: ' ME.message];
            end
            obj.IsDirty = ~info.exists || pythonChanged || gptChanged;
            obj.MapBase = phase_preprocessing_beta.mapBase();
            obj.MapFootprints = phase_preprocessing_beta.collectFootprints(obj.RootDir);
            if info.exists, obj.MapPolygon = bboxPolygon(obj.Config);
            else, obj.MapPolygon = zeros(0, 2); end
            obj.SlaveFiles = phase_preprocessing_beta.scanSlaves(obj.RootDir);
            obj.Downloader = phase_preprocessing_beta.defaultDownloader(obj.RootDir, obj.MapPolygon);
            obj.UpdateContext = phase_preprocessing_beta.sentinelUpdateContext(obj.RootDir);
            obj.UpdateStatus = obj.UpdateContext.message;
            obj.Transfer = defaultTransfer();

            obj.UIFigure = uifigure('Name', 'PHASE · Preprocessing Beta', ...
                'Color', [1 1 1], 'Position', centeredPosition(1500, 920));
            obj.UIFigure.UserData = obj;
            obj.UIFigure.CloseRequestFcn = @(~,~) delete(obj);
            grid = uigridlayout(obj.UIFigure, [1 1]);
            grid.Padding = [0 0 0 0];
            uiPath = fullfile(obj.RootDir, 'PHASE_Preprocessing', ...
                'phase_preprocessing_beta_ui', 'index.html');
            obj.HTML = uihtml(grid, 'HTMLSource', uiPath);
            obj.HTML.Layout.Row = 1; obj.HTML.Layout.Column = 1;
            obj.HTML.HTMLEventReceivedFcn = @(~,event) obj.onHtmlEvent(event);
            drawnow;
            obj.sendState();

            try
                obj.appendLog(pythonNotice);
                obj.appendLog(gptNotice);
                obj.appendLog('Initialising the proven preprocessing engine and map services…');
                obj.Engine = phase_preprocessing_beta.ProcessingEngine();
                obj.Engine.ExternalLogCallback = @(message) obj.appendLog(message);
                obj.Engine.ExternalProgressCallback = @(progress) obj.onEngineProgress(progress);
                obj.Engine.UIFigure.CloseRequestFcn = @(~,~) obj.hideEngine();
                phase_preprocessing_beta.themeLegacyEngine(obj.Engine);
                phase_preprocessing_beta.applyConfig(obj.Engine, obj.Config);
                obj.Status = ternary(obj.IsDirty, 'idle', 'saved');
                obj.StatusDetail = ternary(obj.IsDirty, ...
                    'Review and save the initial configuration', ...
                    'Configuration loaded; engine ready');
                obj.appendLog('Preprocessing engine ready.');
            catch ME
                obj.Status = 'error'; obj.StatusDetail = ME.message;
                obj.appendLog(['Engine initialisation failed [' ME.identifier ']: ' ME.message]);
                obj.showError('PHASE Preprocessing Beta', ME.message);
            end
            obj.sendState();
        end

        function appendLog(obj, message)
            message = char(string(message));
            timestamp = char(datetime('now', 'Format', 'HH:mm:ss'));
            entry = struct('time', timestamp, 'message', message);
            obj.Logs{end+1} = entry;
            if ~isempty(obj.HTML) && isvalid(obj.HTML)
                try
                    sendEventToHTMLSource(obj.HTML, 'PhaseLog', entry);
                    drawnow limitrate
                catch
                end
            end
        end

        function delete(obj)
            try, obj.stopDownloadTransfer(true); catch, end
            try
                if ~isempty(obj.Engine) && isvalid(obj.Engine)
                    obj.Engine.UIFigure.CloseRequestFcn = [];
                    delete(obj.Engine);
                end
            catch
            end
            try
                if ~isempty(obj.UIFigure) && isvalid(obj.UIFigure)
                    obj.UIFigure.CloseRequestFcn = [];
                    obj.UIFigure.UserData = [];
                    delete(obj.UIFigure);
                end
            catch
            end
        end
    end

    methods (Access = private)
        function onHtmlEvent(obj, event)
            try
                name = lower(char(string(event.HTMLEventName)));
                payload = event.HTMLEventData;
                switch name
                    case 'ready', obj.sendState();
                    case 'load', obj.loadFromDisk();
                    case 'save', obj.saveFromPayload(payload);
                    case 'start', obj.startFromPayload(payload);
                    case 'stop', obj.stopProcessing();
                    case 'browse', obj.browseFromPayload(payload);
                    case 'openworkdir', openFolder(fullfile(obj.RootDir, 'PHASE_Preprocessing'));
                    case 'importimages', obj.updateFromPayload(payload); obj.importImages();
                    case 'openslavesfolder', obj.updateFromPayload(payload); obj.openSlavesFolder();
                    case 'refreshslaves', obj.refreshLocalData('Local image inventory refreshed.');
                    case 'mapaoichanged', obj.mapAoiChanged(payload);
                    case 'refreshmap', obj.refreshMap();
                    case 'maptilesrequested', obj.cacheMapTiles(payload);
                    case 'downloadaoichanged', obj.downloadAoiChanged(payload);
                    case 'downloadsearch', obj.searchAsf(payload);
                    case 'downloadlogin', obj.loginAsf(payload);
                    case 'downloadlogout', obj.logoutAsf();
                    case 'downloadselected', obj.downloadSelectedAsf(payload);
                    case 'stopdownload', obj.stopDownloadTransfer(false);
                    case 'refreshupdate', obj.refreshUpdateContext();
                    case 'searchupdate', obj.searchUpdateImages(payload);
                    case 'downloadupdate', obj.downloadUpdateImages(payload);
                    otherwise, obj.appendLog(['Unknown interface event: ' name]);
                end
            catch ME
                obj.IsRunning = false;
                obj.Status = 'error'; obj.StatusDetail = ME.message;
                obj.appendLog(['Interface action failed [' ME.identifier ']: ' ME.message]);
                obj.showError('PHASE Preprocessing Beta', ME.message);
                obj.sendState();
            end
        end

        function loadFromDisk(obj)
            if obj.IsRunning, return; end
            [cfg, info] = phase_preprocessing_beta.loadConfig(obj.RootDir);
            if ~info.exists
                error('PHASE_Preprocessing_beta:configMissing', ...
                    'No PHASE_Preprocessing/input_preprocessing.mat file exists yet.');
            end
            savedCfg = cfg;
            [cfg.python, pythonInfo] = phase_preprocessing_beta.resolvePython(cfg.python);
            [cfg.gptbin_path,gptInfo] = ...
                phase_preprocessing_beta.resolveGpt(cfg.gptbin_path);
            obj.Config = cfg; obj.SavedConfig = savedCfg;
            obj.IsDirty = pythonInfo.changed || gptInfo.changed;
            obj.MapPolygon = bboxPolygon(obj.Config);
            obj.applyToEngine();
            obj.Status = ternary(obj.IsDirty, 'idle', 'saved');
            obj.StatusDetail = ternary(obj.IsDirty, ...
                'Python 3 detected; save the resolved executable path', ...
                'Configuration loaded');
            obj.appendLog(['Loaded configuration: ' info.path]);
            obj.appendLog(sprintf('Python %s selected from %s: %s', ...
                pythonInfo.version, pythonInfo.source, cfg.python));
            if gptInfo.exists
                obj.appendLog(['SNAP gpt selected: ' cfg.gptbin_path]);
            end
            obj.sendState();
        end

        function saveFromPayload(obj, payload)
            if obj.IsRunning, return; end
            candidate = obj.configFromPayload(payload);
            [candidate.python, pythonInfo] = ...
                phase_preprocessing_beta.resolvePython(candidate.python);
            [candidate.gptbin_path,gptInfo] = ...
                phase_preprocessing_beta.resolveGpt(candidate.gptbin_path);
            [errors, warnings] = phase_preprocessing_beta.validateConfig(candidate, false);
            if ~isempty(errors)
                error('PHASE_Preprocessing_beta:invalidConfiguration', '%s', strjoin(errors, newline));
            end
            if bboxChanged(obj.Config, candidate), obj.MapPolygon = bboxPolygon(candidate); end
            pathValue = phase_preprocessing_beta.saveConfig(obj.RootDir, candidate);
            obj.Config = candidate; obj.SavedConfig = candidate; obj.IsDirty = false;
            obj.applyToEngine();
            obj.Status = 'saved'; obj.StatusDetail = 'Configuration saved and ready';
            obj.appendLog(['Saved configuration: ' pathValue]);
            obj.appendLog(sprintf('Python %s verified: %s', ...
                pythonInfo.version, candidate.python));
            if gptInfo.exists
                obj.appendLog(['SNAP gpt verified: ' candidate.gptbin_path]);
            end
            for k = 1:numel(warnings), obj.appendLog(['Warning: ' warnings{k}]); end
            obj.sendState();
        end

        function startFromPayload(obj, payload)
            if obj.IsRunning, return; end
            obj.requireEngine();
            candidate = obj.configFromPayload(payload);
            candidate.python = phase_preprocessing_beta.resolvePython(candidate.python);
            candidate.gptbin_path = ...
                phase_preprocessing_beta.resolveGpt(candidate.gptbin_path);
            [errors, warnings] = phase_preprocessing_beta.validateConfig(candidate, true);
            if ~isempty(errors)
                error('PHASE_Preprocessing_beta:invalidConfiguration', '%s', strjoin(errors, newline));
            end
            if isempty(obj.SavedConfig) || ...
                    ~phase_preprocessing_beta.configsEqual(candidate, obj.SavedConfig)
                error('PHASE_Preprocessing_beta:unsavedConfiguration', ...
                    'The visible configuration differs from input_preprocessing.mat. Press Save before Start.');
            end
            obj.Config = candidate; obj.applyToEngine();
            obj.Logs = {}; obj.IsRunning = true; obj.StopRequested = false;
            obj.resetRunProgress();
            obj.Status = 'running';
            obj.StatusDetail = sprintf('%s preprocessing from step %d', ...
                constellationLabel(candidate), candidate.first_step);
            obj.sendState();
            for k = 1:numel(warnings), obj.appendLog(['Warning: ' warnings{k}]); end
            obj.appendLog(['Starting ' constellationLabel(candidate) ' preprocessing.']);
            previous = pwd; cleanup = onCleanup(@() cd(previous)); %#ok<NASGU>
            cd(obj.RootDir); drawnow;
            mkdirWarning = warning('off','MATLAB:MKDIR:DirectoryExists');
            warningCleanup = onCleanup(@() warning(mkdirWarning)); %#ok<NASGU>
            try
                obj.Engine.StartButtonPushed([]);
                obj.IsRunning = false;
                phase_preprocessing_beta.saveConfig(obj.RootDir,candidate);
                if obj.StopRequested
                    obj.Status = 'idle';
                    obj.StatusDetail = 'Preprocessing stopped; cleanup was skipped';
                    obj.finishRunProgress('Stopped after current step',obj.RunProgress.percentage);
                    obj.appendLog('Preprocessing stopped. Selected cleanup actions were not applied.');
                else
                    obj.finalizeBetaProducts(candidate);
                    obj.Status = 'success';
                    obj.StatusDetail = 'Preprocessing completed';
                    obj.finishRunProgress('Preprocessing completed',100);
                end
            catch ME
                try, phase_preprocessing_beta.saveConfig(obj.RootDir,candidate); catch, end
                obj.IsRunning = false;
                if obj.StopRequested || strcmp(ME.identifier,'PHASE:ProcessingStopped')
                    obj.Status = 'idle';
                    obj.StatusDetail = 'Processing force-stopped; cleanup was skipped';
                    obj.finishRunProgress('Processing force-stopped',obj.RunProgress.percentage);
                    obj.appendLog('Processing was force-stopped. Partial products from the interrupted SNAP operation may be incomplete; cleanup was skipped.');
                else
                    obj.Status = 'error'; obj.StatusDetail = ME.message;
                    obj.finishRunProgress('Processing failed',obj.RunProgress.percentage);
                    obj.appendLog(['Processing failed [' ME.identifier ']: ' ME.message]);
                end
            end
            obj.sendState();
        end

        function stopProcessing(obj)
            if isempty(obj.Engine) || ~isvalid(obj.Engine), return; end
            answer = uiconfirm(obj.UIFigure, ...
                ['Force-stop the active preprocessing operation now?' newline newline ...
                 'MATLAB will terminate the complete CMD/Python/SNAP process tree. ' ...
                 'The product currently being written may be incomplete and should be recreated.'], ...
                'Force stop preprocessing', ...
                'Options',{'Force stop now','Cancel'}, ...
                'DefaultOption','Cancel','CancelOption','Cancel');
            if ~strcmp(answer,'Force stop now'), return; end
            obj.StopRequested = true;
            obj.Engine.StopFlag = true;
            stopped = false;
            try, stopped = obj.Engine.forceStopActiveProcess(); catch, end
            if stopped
                obj.StatusDetail = 'Force stop sent to CMD, Python and SNAP';
                obj.appendLog('Force stop sent: the active CMD/Python/SNAP process tree was terminated.');
            else
                obj.StatusDetail = 'Force stop requested; interrupting at the next engine yield';
                obj.appendLog('Force stop requested. No active external process was found; PHASE will stop at the next engine yield.');
            end
            obj.sendState();
        end

        function browseFromPayload(obj, payload)
            if obj.IsRunning, return; end
            obj.updateFromPayload(payload);
            if ~isstruct(payload) || ~isfield(payload, 'field')
                error('PHASE_Preprocessing_beta:missingBrowseField', 'Browse action has no target field.');
            end
            fieldName = char(string(payload.field));
            if any(strcmp(fieldName, {'dem_file','dem_file_coreg'}))
                [file, folder] = uigetfile({'*.tif;*.tiff','GeoTIFF DEM (*.tif, *.tiff)'}, 'Select external DEM');
                if isequal(file, 0), return; end
                selected = fullfile(folder, file);
            elseif strcmp(fieldName, 'gptbin_path')
                [file, folder] = uigetfile({'gpt*','SNAP gpt executable'; '*.*','All files'}, 'Select SNAP gpt executable');
                if isequal(file, 0), return; end
                selected = fullfile(folder, file);
            else
                error('PHASE_Preprocessing_beta:invalidBrowseField', 'Unsupported path field: %s', fieldName);
            end
            obj.Config.(fieldName) = selected; obj.IsDirty = true;
            obj.Status = 'idle'; obj.StatusDetail = 'Unsaved changes'; obj.sendState();
        end

        function importImages(obj)
            if obj.IsRunning, return; end
            isSEN = strcmp(obj.Config.constellation, 'SEN');
            if isSEN
                [files, sourceFolder] = uigetfile({'*.zip','Sentinel-1 ZIP files (*.zip)'}, ...
                    'Select Sentinel-1 images', 'MultiSelect', 'on');
            else
                [files, sourceFolder] = uigetfile({'*.h5','COSMO-SkyMed HDF5 files (*.h5)'}, ...
                    'Select COSMO-SkyMed images', 'MultiSelect', 'on');
            end
            if isequal(files,0), return; end
            if ischar(files), files = {files}; end

            mode = 'Copy';
            if isSEN
                mode = uiconfirm(obj.UIFigure, ...
                    ['Move or copy the selected Sentinel-1 ZIP files into PHASE_Preprocessing/slaves?' newline newline ...
                    'Move avoids duplicating large files. Copy keeps the originals.'], ...
                    'Import Sentinel-1 images', 'Options',{'Move','Copy','Cancel'}, ...
                    'DefaultOption','Move','CancelOption','Cancel');
                if strcmp(mode,'Cancel'), return; end
            end

            destination = obj.slavesFolder();
            if ~isfolder(destination), mkdir(destination); end
            dialog = uiprogressdlg(obj.UIFigure,'Title','Importing images', ...
                'Message','Preparing import…','Value',0);
            cleanup = onCleanup(@() closeProgress(dialog)); %#ok<NASGU>
            imported = 0; skipped = 0; failures = {};
            for k = 1:numel(files)
                dialog.Value = k/numel(files);
                dialog.Message = sprintf('%s %d of %d: %s',mode,k,numel(files),files{k});
                source = fullfile(sourceFolder,files{k});
                target = fullfile(destination,files{k});
                if isfile(target), skipped = skipped + 1; continue; end
                try
                    if strcmp(mode,'Move'), [ok,message] = movefile(source,target);
                    else, [ok,message] = copyfile(source,target); end
                    if ~ok, error('PHASE:ImageImportFailed','%s',message); end
                    imported = imported + 1;
                catch ME
                    failures{end+1} = sprintf('%s: %s',files{k},ME.message); %#ok<AGROW>
                end
            end
            obj.refreshLocalData(sprintf('Imported %d image(s); skipped %d existing file(s).',imported,skipped));
            if ~isempty(failures), obj.showError('Import failed',strjoin(failures,newline)); end
        end

        function openSlavesFolder(obj)
            folder = obj.slavesFolder();
            if ~isfolder(folder), mkdir(folder); end
            openFolder(folder);
        end

        function refreshLocalData(obj, message)
            obj.SlaveFiles = phase_preprocessing_beta.scanSlaves(obj.RootDir);
            obj.MapFootprints = phase_preprocessing_beta.collectFootprints(obj.RootDir);
            obj.UpdateContext = phase_preprocessing_beta.sentinelUpdateContext(obj.RootDir);
            if nargin > 1 && ~isempty(message), obj.appendLog(message); end
            obj.sendState();
        end

        function folder = slavesFolder(obj)
            folder = fullfile(obj.RootDir,'PHASE_Preprocessing','slaves');
        end

        function downloadAoiChanged(obj, payload)
            if ~isstruct(payload) || ~isfield(payload,'polygon')
                error('PHASE:InvalidDownloaderAOI','The downloader map did not provide a polygon.');
            end
            polygon = numericMatrix(payload.polygon);
            if size(polygon,1) < 3 || size(polygon,2) ~= 2
                error('PHASE:InvalidDownloaderAOI','Draw at least three valid AOI vertices.');
            end
            if ~isequal(polygon(1,:),polygon(end,:)), polygon(end+1,:) = polygon(1,:); end
            obj.Downloader.polygon = polygon;
            obj.Downloader.status = sprintf('Download AOI ready with %d vertices.',size(polygon,1)-1);
            obj.sendState();
        end

        function searchAsf(obj, payload)
            if isstruct(payload) && isfield(payload,'config'), obj.updateFromPayload(payload); end
            if ~isstruct(payload) || ~isfield(payload,'filters')
                error('PHASE:MissingASFFilters','Downloader filters were not provided.');
            end
            if size(obj.Downloader.polygon,1) < 4
                error('PHASE:MissingDownloaderAOI','Draw a download AOI before searching ASF.');
            end
            filters = normalizeFilters(payload.filters,obj.Downloader.filters);
            request = buildSearchRequest(filters,obj.Downloader.polygon);
            obj.Downloader.filters = filters;
            obj.Downloader.busy = true; obj.Downloader.progress = 0;
            obj.Downloader.status = 'Searching ASF…'; obj.sendState();
            try
                asfFolder = fullfile(obj.RootDir,'downloadasf');
                deleteIfExists(fullfile(asfFolder,'search_summary.json'));
                deleteIfExists(fullfile(asfFolder,'download_data.json'));
                writeJson(fullfile(asfFolder,'search_request.json'),request);
                controller = fullfile(obj.RootDir,'downloadasf','controller.py');
                [status,output] = obj.runPython(controller,{'search'});
                if status ~= 0, error('PHASE:ASFSearchFailed','%s',strtrim(output)); end
                if ~isfile(fullfile(asfFolder,'search_summary.json')) || ...
                        ~isfile(fullfile(asfFolder,'download_data.json'))
                    error('PHASE:ASFSearchFailed','The ASF backend returned no search result files. %s',strtrim(output));
                end
                backupAsfFiles(obj.RootDir);
                search = phase_preprocessing_beta.readAsfSearch(obj.RootDir);
                obj.Downloader.results = search.results;
                obj.Downloader.count = search.count;
                obj.Downloader.totalSizeGB = search.totalSizeGB;
                obj.Downloader.recommended = search.recommended;
                obj.Downloader.status = sprintf('ASF returned %d product(s), %.2f GB.', ...
                    search.count,search.totalSizeGB);
                obj.Downloader.busy = false; obj.Downloader.progress = 100;
                obj.appendLog(obj.Downloader.status); obj.sendState();
            catch ME
                obj.Downloader.busy = false; obj.Downloader.status = ['Search failed: ' ME.message];
                obj.sendState(); rethrow(ME);
            end
        end

        function loginAsf(obj, payload)
            if ~isstruct(payload) || ~isfield(payload,'username') || ~isfield(payload,'password')
                error('PHASE:MissingEarthdataCredentials','Earthdata username and password are required.');
            end
            username = strtrim(char(string(payload.username)));
            password = char(string(payload.password));
            if isempty(username) || isempty(password)
                error('PHASE:MissingEarthdataCredentials','Earthdata username and password are required.');
            end
            obj.Downloader.busy = true; obj.Downloader.status = 'Checking Earthdata credentials…'; obj.sendState();
            try
                writeJson(fullfile(obj.RootDir,'downloadasf','login_request.json'), ...
                    struct('username',username,'password',password));
                controller = fullfile(obj.RootDir,'downloadasf','controller.py');
                [status,output] = obj.runPython(controller,{'login'});
                resultPath = fullfile(obj.RootDir,'downloadasf','login_result.json');
                if status ~= 0 || ~isfile(resultPath)
                    error('PHASE:EarthdataLoginFailed','%s',strtrim(output));
                end
                result = jsondecode(fileread(resultPath));
                if ~isfield(result,'status') || ~strcmp(string(result.status),"success")
                    error('PHASE:EarthdataLoginFailed','Earthdata rejected the credentials.');
                end
                obj.Downloader.loggedIn = true; obj.Downloader.username = username;
                obj.Downloader.status = ['Signed in as ' username '.']; obj.Downloader.busy = false;
                obj.appendLog(['Earthdata login verified for ' username '.']); obj.sendState();
            catch ME
                obj.Downloader.loggedIn = false; obj.Downloader.username = '';
                obj.Downloader.busy = false; obj.Downloader.status = ['Login failed: ' ME.message];
                obj.sendState(); rethrow(ME);
            end
        end

        function logoutAsf(obj)
            deleteIfExists(fullfile(obj.RootDir,'downloadasf','login_result.json'));
            deleteIfExists(fullfile(obj.RootDir,'downloadasf','login_request.json'));
            obj.Downloader.loggedIn = false; obj.Downloader.username = '';
            obj.Downloader.status = 'Signed out from Earthdata.'; obj.sendState();
        end

        function downloadSelectedAsf(obj, payload)
            if ~obj.Downloader.loggedIn
                error('PHASE:EarthdataLoginRequired','Sign in to Earthdata before downloading.');
            end
            names = payloadNames(payload);
            products = obj.Downloader.results;
            if isempty(products) || isempty(names)
                error('PHASE:NoASFSelection','Select at least one ASF result.');
            end
            selected = products(ismember(string({products.sceneName}),names));
            validateCompatibleSelection(selected);
            existing = countExisting(obj.slavesFolder(),selected);
            message = sprintf(['Download %d selected Sentinel-1 image(s)?\n\nTotal size: %.2f GB\n' ...
                'Already in slaves: %d\nNew downloads: %d\n\n' ...
                'Completed files are kept, interrupted .part files can be resumed, ' ...
                'and no other ZIP in slaves will be removed.'], ...
                numel(selected),sum([selected.sizeGB]),existing,numel(selected)-existing);
            answer = uiconfirm(obj.UIFigure,message,'Confirm ASF download', ...
                'Options',{'Download','Cancel'},'DefaultOption','Download','CancelOption','Cancel');
            if ~strcmp(answer,'Download'), return; end

            saveAsfSelection(obj.RootDir,selected);
            entries = initialDownloadEntries(obj.RootDir, selected);
            obj.startDownloadTransfer('initial', entries);
        end

        function refreshUpdateContext(obj)
            obj.UpdateContext = phase_preprocessing_beta.sentinelUpdateContext(obj.RootDir);
            obj.UpdateStatus = obj.UpdateContext.message;
            obj.sendState();
        end

        function searchUpdateImages(obj, payload)
            obj.UpdateContext = phase_preprocessing_beta.sentinelUpdateContext(obj.RootDir);
            if ~obj.UpdateContext.available, error('PHASE:UpdateContextMissing','%s',obj.UpdateContext.message); end
            if ~isstruct(payload) || ~isfield(payload,'endDate')
                error('PHASE:UpdateEndDateMissing','Choose the update search end date.');
            end
            endDate = regexprep(char(string(payload.endDate)),'[^0-9]','');
            if numel(endDate) ~= 8, error('PHASE:InvalidUpdateEndDate','Use a valid update end date.'); end
            if datetime(endDate,'InputFormat','yyyyMMdd') < datetime(obj.UpdateContext.latestDate,'InputFormat','yyyyMMdd')
                error('PHASE:InvalidUpdateRange','The end date cannot precede the latest local acquisition.');
            end
            obj.UpdateBusy = true; obj.UpdateStatus = 'Searching ASF for newer compatible scenes…'; obj.sendState();
            outputPath = [tempname '.json']; cleanup = onCleanup(@() deleteIfExists(outputPath)); %#ok<NASGU>
            script = fullfile(obj.RootDir,'pythonScripts','search_update_sentinel1_images.py');
            arguments = {'--reference',obj.UpdateContext.referenceZip, ...
                '--start-date',obj.UpdateContext.latestDate,'--end-date',endDate, ...
                '--exclude-dates',strjoin(obj.UpdateContext.localDates,','), ...
                '--output',outputPath};
            try
                [status,output] = obj.runPython(script,arguments);
                if status ~= 0 || ~isfile(outputPath), error('PHASE:UpdateSearchFailed','%s',strtrim(output)); end
                data = jsondecode(fileread(outputPath));
                if isfield(data,'results'), obj.UpdateResults = normalizeUpdateResults(data.results);
                else, obj.UpdateResults = struct([]); end
                obj.UpdateBusy = false;
                obj.UpdateStatus = sprintf('Found %d newer compatible image(s).',numel(obj.UpdateResults));
                obj.appendLog(obj.UpdateStatus); obj.sendState();
            catch ME
                obj.UpdateBusy = false; obj.UpdateStatus = ['Update search failed: ' ME.message];
                obj.sendState(); rethrow(ME);
            end
        end

        function downloadUpdateImages(obj, payload)
            if ~obj.Downloader.loggedIn
                error('PHASE:EarthdataLoginRequired','Sign in to Earthdata before downloading update images.');
            end
            names = payloadNames(payload);
            if isempty(names) || isempty(obj.UpdateResults)
                error('PHASE:NoUpdateSelection','Select at least one update image.');
            end
            selected = obj.UpdateResults(ismember(string({obj.UpdateResults.sceneName}),names));
            if isempty(selected), error('PHASE:NoUpdateSelection','Select at least one update image.'); end
            entries = updateDownloadEntries(selected);
            obj.startDownloadTransfer('update', entries);
        end

        function startDownloadTransfer(obj, kind, entries)
            if isstruct(obj.Transfer) && isfield(obj.Transfer,'active') && obj.Transfer.active
                error('PHASE:DownloadAlreadyRunning', ...
                    'A PHASE image download is already running. Stop it before starting another.');
            end
            if isempty(entries)
                error('PHASE:NoDownloadEntries','The selected images contain no download URLs.');
            end

            pythonExecutable = phase_preprocessing_beta.resolvePython(obj.Config.python);
            obj.TransferDirectory = tempname;
            mkdir(obj.TransferDirectory);
            manifestPath = fullfile(obj.TransferDirectory,'manifest.json');
            progressPath = fullfile(obj.TransferDirectory,'progress.json');
            resultPath = fullfile(obj.TransferDirectory,'result.json');
            stopPath = fullfile(obj.TransferDirectory,'stop.request');
            obj.TransferLogPath = fullfile(obj.TransferDirectory,'download.log');
            writeJson(manifestPath,struct('files',{entries}));

            script = fullfile(obj.RootDir,'pythonScripts','phase_download_manager.py');
            credentials = fullfile(obj.RootDir,'downloadasf','login_request.json');
            parts = {pythonExecutable,script,'--manifest',manifestPath, ...
                '--destination',obj.slavesFolder(),'--credentials',credentials, ...
                '--progress',progressPath,'--output',resultPath, ...
                '--stop',stopPath,'--kind',kind};
            command = java.util.ArrayList();
            for k = 1:numel(parts)
                command.add(java.lang.String(char(string(parts{k}))));
            end

            try
                builder = java.lang.ProcessBuilder(command);
                builder.directory(java.io.File(obj.RootDir));
                builder.redirectErrorStream(true);
                builder.redirectOutput(java.io.File(obj.TransferLogPath));
                obj.TransferProcess = builder.start();
            catch ME
                cleanupTransferDirectory(obj.TransferDirectory);
                obj.TransferDirectory = '';
                error('PHASE:DownloadLaunchFailed', ...
                    'Could not launch the background downloader: %s',ME.message);
            end

            obj.Transfer = defaultTransfer();
            obj.Transfer.active = true;
            obj.Transfer.canStop = true;
            obj.Transfer.kind = kind;
            obj.Transfer.phase = 'preparing';
            obj.Transfer.totalFiles = numel(entries);
            obj.Transfer.message = sprintf('Preparing %d image download(s)…',numel(entries));
            if strcmp(kind,'initial')
                obj.Downloader.busy = true;
                obj.Downloader.progress = 0;
                obj.Downloader.status = obj.Transfer.message;
            else
                obj.UpdateBusy = true;
                obj.UpdateStatus = obj.Transfer.message;
            end

            obj.TransferTimer = timer('ExecutionMode','fixedSpacing', ...
                'Period',0.4,'BusyMode','drop', ...
                'TimerFcn',@(~,~) obj.pollDownloadTransfer());
            start(obj.TransferTimer);
            obj.appendLog(sprintf('Started background %s download of %d image(s).', ...
                kind,numel(entries)));
            obj.sendState();
        end

        function pollDownloadTransfer(obj)
            if ~isstruct(obj.Transfer) || ~obj.Transfer.active, return; end
            previousCompleted = obj.Transfer.completedFiles;
            progressPath = fullfile(obj.TransferDirectory,'progress.json');
            if isfile(progressPath)
                try
                    progress = jsondecode(fileread(progressPath));
                    obj.Transfer = mergeTransferProgress(obj.Transfer,progress);
                catch
                end
            end
            if obj.Transfer.completedFiles ~= previousCompleted
                obj.SlaveFiles = phase_preprocessing_beta.scanSlaves(obj.RootDir);
                obj.MapFootprints = phase_preprocessing_beta.collectFootprints(obj.RootDir);
                obj.UpdateContext = phase_preprocessing_beta.sentinelUpdateContext(obj.RootDir);
            end
            if strcmp(obj.Transfer.kind,'initial')
                obj.Downloader.progress = obj.Transfer.percentage;
                obj.Downloader.status = obj.Transfer.message;
            else
                obj.UpdateStatus = obj.Transfer.message;
            end

            alive = false;
            try, alive = obj.TransferProcess.isAlive(); catch, end
            if alive
                obj.sendState();
            else
                obj.finishDownloadTransfer();
            end
        end

        function finishDownloadTransfer(obj)
            resultPath = fullfile(obj.TransferDirectory,'result.json');
            result = struct('status','failed','failedCount',1);
            if isfile(resultPath)
                try, result = jsondecode(fileread(resultPath)); catch, end
            end
            status = char(string(result.status));
            countFields = {'completedFiles','downloadedCount','skippedCount','failedCount'};
            targetFields = {'completedFiles','downloadedFiles','skippedFiles','failedFiles'};
            for k = 1:numel(countFields)
                if isfield(result,countFields{k})
                    obj.Transfer.(targetFields{k}) = double(result.(countFields{k}));
                end
            end
            obj.Transfer.active = false;
            obj.Transfer.canStop = false;
            obj.Transfer.phase = status;
            if strcmp(status,'completed')
                obj.Transfer.percentage = 100;
                obj.Transfer.message = sprintf('Completed %d of %d images.', ...
                    obj.Transfer.completedFiles,obj.Transfer.totalFiles);
            elseif strcmp(status,'stopped')
                obj.Transfer.message = sprintf('Stopped after %d of %d images.', ...
                    obj.Transfer.completedFiles,obj.Transfer.totalFiles);
            else
                failedCount = 1;
                if isfield(result,'failedCount'), failedCount = double(result.failedCount); end
                obj.Transfer.failedFiles = max(obj.Transfer.failedFiles,failedCount);
                obj.Transfer.message = sprintf('%d image download(s) failed.',failedCount);
                if isfile(obj.TransferLogPath)
                    try
                        logText = strtrim(fileread(obj.TransferLogPath));
                        if ~isempty(logText), obj.appendLog(['Downloader: ' logText]); end
                    catch
                    end
                end
            end

            if strcmp(obj.Transfer.kind,'initial')
                obj.Downloader.busy = false;
                obj.Downloader.progress = obj.Transfer.percentage;
                obj.Downloader.status = obj.Transfer.message;
            else
                obj.UpdateBusy = false;
                obj.UpdateStatus = obj.Transfer.message;
            end
            obj.disposeTransferTimer();
            obj.TransferProcess = [];
            obj.refreshLocalData(obj.Transfer.message);
            cleanupTransferDirectory(obj.TransferDirectory);
            obj.TransferDirectory = '';
            obj.TransferLogPath = '';
        end

        function stopDownloadTransfer(obj, closing)
            if nargin < 2, closing = false; end
            if isempty(obj.Transfer) || ~isstruct(obj.Transfer) || ~obj.Transfer.active
                obj.disposeTransferTimer();
                return;
            end
            stopPath = fullfile(obj.TransferDirectory,'stop.request');
            try
                writeJson(stopPath,struct('requested',true));
            catch
            end
            try
                if obj.TransferProcess.isAlive()
                    obj.TransferProcess.destroyForcibly();
                end
            catch
            end
            obj.Transfer.active = false;
            obj.Transfer.canStop = false;
            obj.Transfer.phase = 'stopped';
            obj.Transfer.message = sprintf('Stopped after %d of %d images; partial files were kept.', ...
                obj.Transfer.completedFiles,obj.Transfer.totalFiles);
            if strcmp(obj.Transfer.kind,'initial')
                obj.Downloader.busy = false;
                obj.Downloader.status = obj.Transfer.message;
            else
                obj.UpdateBusy = false;
                obj.UpdateStatus = obj.Transfer.message;
            end
            obj.disposeTransferTimer();
            obj.TransferProcess = [];
            if ~closing
                obj.appendLog(obj.Transfer.message);
                obj.refreshLocalData(obj.Transfer.message);
            end
            cleanupTransferDirectory(obj.TransferDirectory);
            obj.TransferDirectory = '';
            obj.TransferLogPath = '';
        end

        function disposeTransferTimer(obj)
            try
                if ~isempty(obj.TransferTimer) && isvalid(obj.TransferTimer)
                    stop(obj.TransferTimer);
                    delete(obj.TransferTimer);
                end
            catch
            end
            obj.TransferTimer = [];
        end

        function mapAoiChanged(obj, payload)
            if obj.IsRunning, return; end
            if ~isstruct(payload) || ~isfield(payload, 'polygon') || ~isfield(payload, 'bbox')
                error('PHASE_Preprocessing_beta:invalidMapAOI', 'The map did not provide polygon and bounding-box data.');
            end
            polygon = numericMatrix(payload.polygon);
            bbox = payload.bbox;
            required = {'minLon','maxLon','minLat','maxLat'};
            if size(polygon, 1) < 3 || size(polygon, 2) ~= 2 || ...
                    ~all(cellfun(@(name) isfield(bbox, name), required))
                error('PHASE_Preprocessing_beta:invalidMapAOI', 'Draw a polygon with at least three valid vertices.');
            end
            values = cellfun(@(name) double(bbox.(name)), required);
            if any(~isfinite(values)) || values(1) >= values(2) || values(3) >= values(4)
                error('PHASE_Preprocessing_beta:invalidMapBounds', 'The polygon bounding box is invalid.');
            end
            if ~isequal(polygon(1,:), polygon(end,:)), polygon(end+1,:) = polygon(1,:); end
            obj.MapPolygon = polygon;
            obj.Config.lon_min = values(1); obj.Config.lon_max = values(2);
            obj.Config.lat_min = values(3); obj.Config.lat_max = values(4);
            if obj.Config.auto_epsg
                obj.Config.epsg_code = phase_preprocessing_beta.estimateEpsg( ...
                    values(1),values(3),values(2),values(4));
            end
            obj.IsDirty = isempty(obj.SavedConfig) || ...
                ~phase_preprocessing_beta.configsEqual(obj.Config, obj.SavedConfig);
            obj.Status = 'idle'; obj.StatusDetail = 'AOI polygon updated; save the new bounding box';
            if ~isempty(obj.Engine) && isvalid(obj.Engine)
                phase_preprocessing_beta.applyConfig(obj.Engine, obj.Config);
                obj.applyDownloaderAoi();
            end
            obj.appendLog(sprintf('AOI polygon updated: lon %.6f..%.6f, lat %.6f..%.6f.', ...
                values(1), values(2), values(3), values(4)));
            if obj.Config.auto_epsg
                obj.appendLog(sprintf('Output CRS automatically updated to EPSG:%d.', ...
                    obj.Config.epsg_code));
            end
            obj.sendState();
        end

        function refreshMap(obj)
            obj.MapFootprints = phase_preprocessing_beta.collectFootprints(obj.RootDir);
            obj.appendLog(sprintf('Map refreshed with %d footprint(s).', numel(obj.MapFootprints)));
            obj.sendState();
        end

        function cacheMapTiles(obj, payload)
            if obj.MapTileBusy || ~isstruct(payload) || ~isfield(payload,'requests')
                return;
            end
            requests = payload.requests;
            if isempty(requests), return; end
            if numel(requests) > 120
                obj.appendLog('Map tile request ignored because it exceeded the safe batch limit.');
                return;
            end

            obj.MapTileBusy = true;
            requestPath = [tempname '.json']; outputPath = [tempname '.json'];
            cleanup = onCleanup(@() cleanupFiles({requestPath,outputPath})); %#ok<NASGU>
            try
                writeJson(requestPath,struct('requests',{requests}));
                cacheRoot = fullfile(obj.RootDir,'PHASE_Preprocessing', ...
                    'phase_preprocessing_beta_ui','map_tiles');
                if ~isfolder(cacheRoot), mkdir(cacheRoot); end
                script = fullfile(obj.RootDir,'pythonScripts','cache_phase_map_tiles.py');
                [status,output] = obj.runPython(script,{'--cache-root',cacheRoot, ...
                    '--requests',requestPath,'--output',outputPath});
                if status ~= 0 || ~isfile(outputPath)
                    error('PHASE:MapTileCacheFailed','%s',strtrim(output));
                end
                result = jsondecode(fileread(outputPath));
                keys = successfulTileKeys(result);
                failedCount = 0;
                if isfield(result,'failed'), failedCount = numel(result.failed); end
                sendEventToHTMLSource(obj.HTML,'MapTilesReady', ...
                    struct('keys',{keys},'failedCount',failedCount));
                if failedCount > 0
                    obj.appendLog(sprintf('Satellite map cache loaded %d tile(s); %d failed.', ...
                        numel(keys),failedCount));
                end
            catch ME
                obj.appendLog(['Satellite background unavailable; using vector fallback: ' ME.message]);
                try
                    sendEventToHTMLSource(obj.HTML,'MapTilesReady', ...
                        struct('keys',{{}},'failedCount',numel(requests)));
                catch
                end
            end
            obj.MapTileBusy = false;
        end

        function applyDownloaderAoi(obj)
            if isempty(obj.Engine) || ~isvalid(obj.Engine) || size(obj.MapPolygon, 1) < 4
                return;
            end
            obj.Engine.DownloaderAOIType = "polygon";
            obj.Engine.DownloaderPolygonCoords = obj.MapPolygon;
            obj.Engine.DownloaderAOICorners = [];
            obj.Engine.MinLongitudeEditField.Value = obj.Config.lon_min;
            obj.Engine.MaxLongitudeEditField.Value = obj.Config.lon_max;
            obj.Engine.MinLatitudeEditField.Value = obj.Config.lat_min;
            obj.Engine.MaxLatitudeEditField.Value = obj.Config.lat_max;
        end

        function [status, output] = runPython(obj, scriptPath, arguments)
            pythonExecutable = phase_preprocessing_beta.resolvePython(obj.Config.python);
            parts = [{pythonExecutable},{scriptPath},arguments];
            quoted = cellfun(@quoteCommandArgument,parts,'UniformOutput',false);
            command = strjoin(quoted,' ');
            [status,output] = system(command);
        end

        function hideEngine(obj)
            try, obj.Engine.UIFigure.Visible = 'off'; catch, end
        end

        function applyToEngine(obj)
            obj.requireEngine();
            phase_preprocessing_beta.applyConfig(obj.Engine, obj.Config);
        end

        function requireEngine(obj)
            if isempty(obj.Engine) || ~isvalid(obj.Engine) || ...
                    isempty(obj.Engine.UIFigure) || ~isvalid(obj.Engine.UIFigure)
                error('PHASE_Preprocessing_beta:engineUnavailable', ...
                    'The preprocessing engine is not available. Restart the beta app.');
            end
        end

        function updateFromPayload(obj, payload)
            candidate = obj.configFromPayload(payload);
            if bboxChanged(obj.Config, candidate), obj.MapPolygon = bboxPolygon(candidate); end
            obj.Config = candidate;
            obj.IsDirty = isempty(obj.SavedConfig) || ...
                ~phase_preprocessing_beta.configsEqual(obj.Config, obj.SavedConfig);
        end

        function cfg = configFromPayload(obj, payload)
            if ~isstruct(payload) || ~isfield(payload, 'config')
                error('PHASE_Preprocessing_beta:missingConfig', 'The interface did not send a configuration.');
            end
            cfg = phase_preprocessing_beta.uiToConfig(payload.config, obj.Config);
        end

        function finalizeBetaProducts(obj, cfg)
            if cfg.remove_slaves_after_processing
                obj.removeSourceImages();
            end
            if ~cfg.generate_coherence
                obj.removeProcessingFolder('coherence','coherence output disabled');
            end
            if ~cfg.generate_lia
                obj.removeProcessingFolder('lia','LIA output disabled');
            end
            if cfg.remove_split_after_processing
                if strcmp(cfg.constellation,'SEN'), folder = 'split'; else, folder = 'subset'; end
                obj.removeProcessingFolder(folder,'selected SNAP cleanup');
            end
            if cfg.remove_coreg_after_processing
                obj.removeProcessingFolder('coreg','selected SNAP cleanup');
            end
            if cfg.remove_ifg_after_processing
                obj.removeProcessingFolder('ifg','selected SNAP cleanup');
            end
        end

        function removeSourceImages(obj)
            pathValue = fullfile(obj.RootDir,'PHASE_Preprocessing','slaves');
            try
                if isfolder(pathValue), rmdir(pathValue,'s'); end
                mkdir(pathValue);
                obj.appendLog(sprintf('Removed downloaded source images from %s.',pathValue));
            catch ME
                obj.appendLog(sprintf('Warning: could not clean source images in %s: %s', ...
                    pathValue,ME.message));
            end
        end

        function removeProcessingFolder(obj, name, reason)
            pathValue = fullfile(obj.RootDir,'PHASE_Preprocessing',name);
            if ~isfolder(pathValue), return; end
            try
                rmdir(pathValue,'s');
                obj.appendLog(sprintf('Removed %s (%s).',pathValue,reason));
            catch ME
                obj.appendLog(sprintf('Warning: could not remove %s: %s',pathValue,ME.message));
            end
        end

        function resetRunProgress(obj)
            obj.RunStartedAt = datetime('now');
            obj.RunProgress = struct( ...
                'percentage',0, ...
                'phase','Preparing preprocessing', ...
                'elapsedSeconds',0, ...
                'etaSeconds',-1, ...
                'startedAt',char(datetime(obj.RunStartedAt,'Format',"yyyy-MM-dd'T'HH:mm:ss")), ...
                'indeterminate',false);
            obj.emitRunProgress();
        end

        function onEngineProgress(obj, progress)
            if ~obj.IsRunning || isempty(obj.RunStartedAt), return; end
            percentage = obj.RunProgress.percentage;
            phase = obj.RunProgress.phase;
            indeterminate = false;
            if isstruct(progress)
                if isfield(progress,'percentage')
                    percentage = max(percentage,double(progress.percentage));
                end
                if isfield(progress,'phase')
                    phase = char(string(progress.phase));
                end
                if isfield(progress,'indeterminate')
                    indeterminate = logical(progress.indeterminate);
                end
            end
            elapsed = max(0,seconds(datetime('now')-obj.RunStartedAt));
            eta = -1;
            if percentage >= 8 && percentage < 99 && elapsed >= 15
                eta = elapsed*(100-percentage)/percentage;
            elseif percentage >= 99
                eta = 0;
            end
            obj.RunProgress.percentage = max(0,min(100,percentage));
            obj.RunProgress.phase = phase;
            obj.RunProgress.elapsedSeconds = elapsed;
            obj.RunProgress.etaSeconds = eta;
            obj.RunProgress.indeterminate = indeterminate;
            obj.StatusDetail = phase;
            obj.emitRunProgress();
        end

        function finishRunProgress(obj, phase, percentage)
            if isempty(obj.RunStartedAt), return; end
            elapsed = max(0,seconds(datetime('now')-obj.RunStartedAt));
            obj.RunProgress.percentage = max(0,min(100,double(percentage)));
            obj.RunProgress.phase = char(string(phase));
            obj.RunProgress.elapsedSeconds = elapsed;
            if obj.RunProgress.percentage >= 100
                obj.RunProgress.etaSeconds = 0;
            else
                obj.RunProgress.etaSeconds = -1;
            end
            obj.RunProgress.indeterminate = false;
            obj.emitRunProgress();
        end

        function emitRunProgress(obj)
            if isempty(obj.HTML) || ~isvalid(obj.HTML), return; end
            try
                sendEventToHTMLSource(obj.HTML,'PhaseRunProgress',obj.RunProgress);
                drawnow limitrate
            catch
            end
        end

        function sendState(obj)
            if isempty(obj.HTML) || ~isvalid(obj.HTML), return; end
            mapState = struct('coastlines',{obj.MapBase}, ...
                'footprints',{obj.MapFootprints},'polygon',obj.MapPolygon);
            updateState = struct('context',obj.UpdateContext, ...
                'results',{obj.UpdateResults},'status',obj.UpdateStatus,'busy',obj.UpdateBusy);
            state = struct('kind','state','version','0.1.0-beta', ...
                'workDir',fullfile(obj.RootDir,'PHASE_Preprocessing'), ...
                'configPath',fullfile(obj.RootDir,'PHASE_Preprocessing','input_preprocessing.mat'), ...
                'schema',phase_preprocessing_beta.schema(), ...
                'config',phase_preprocessing_beta.configToUi(obj.Config), ...
                'detectedFields',{{}},'dirty',obj.IsDirty,'running',obj.IsRunning, ...
                'status',obj.Status,'statusDetail',obj.StatusDetail,'logs',{obj.Logs}, ...
                'runProgress',obj.RunProgress, ...
                'map',mapState,'slaves',{obj.SlaveFiles}, ...
                'downloader',obj.Downloader,'update',updateState, ...
                'transfer',obj.Transfer);
            obj.HTML.Data = state; drawnow limitrate
        end

        function showError(obj, titleText, message)
            try, uialert(obj.UIFigure, char(string(message)), titleText);
            catch, warning('%s: %s', titleText, char(string(message))); end
        end
    end
end

function pos = centeredPosition(width, height)
screen = get(groot, 'ScreenSize');
width = min(width, max(1050, screen(3)-60)); height = min(height, max(720, screen(4)-100));
pos = [max(20,(screen(3)-width)/2), max(40,(screen(4)-height)/2), width, height];
end

function label = constellationLabel(cfg)
if strcmp(cfg.constellation, 'CSK'), label = 'COSMO-SkyMed'; else, label = 'Sentinel-1'; end
end

function out = ternary(condition, yesValue, noValue)
if condition, out = yesValue; else, out = noValue; end
end

function state = defaultTransfer()
state = struct( ...
    'active',false,'canStop',false,'kind','','phase','idle', ...
    'percentage',0,'currentIndex',0,'totalFiles',0,'currentFile','', ...
    'currentBytes',0,'currentTotalBytes',0,'completedBytes',0, ...
    'expectedTotalBytes',0,'completedFiles',0,'downloadedFiles',0, ...
    'skippedFiles',0,'failedFiles',0,'message','No image download is running.');
end

function state = mergeTransferProgress(state, progress)
fields = {'kind','phase','percentage','currentIndex','totalFiles', ...
    'currentFile','currentBytes','currentTotalBytes','completedBytes', ...
    'expectedTotalBytes','completedFiles','skippedFiles','failedFiles','message'};
for k = 1:numel(fields)
    name = fields{k};
    if isfield(progress,name), state.(name) = progress.(name); end
end
state.canStop = true;
end

function entries = initialDownloadEntries(rootDir, selected)
pathValue = fullfile(rootDir,'downloadasf','download_data.json');
if ~isfile(pathValue)
    error('PHASE:MissingASFDownloadData', ...
        'The selected ASF download metadata file was not created.');
end
data = jsondecode(fileread(pathValue));
if ~isfield(data,'information') || isempty(data.information)
    error('PHASE:MissingASFDownloadData','The ASF download metadata contains no URLs.');
end
information = data.information;
template = struct('name','','url','','sizeBytes',0);
entries = repmat(template,0,1);
for k = 1:numel(selected)
    scene = char(string(selected(k).sceneName));
    match = find(strcmp(string({information.sceneName}),string(scene)),1);
    if isempty(match) || ~isfield(information(match),'url'), continue; end
    entry = template;
    entry.name = [scene '.zip'];
    entry.url = char(string(information(match).url));
    entry.sizeBytes = double(selected(k).sizeBytes);
    entries(end+1,1) = entry; %#ok<AGROW>
end
if numel(entries) ~= numel(selected)
    error('PHASE:MissingASFDownloadURL', ...
        'One or more selected ASF products have no download URL.');
end
end

function entries = updateDownloadEntries(selected)
template = struct('name','','url','','sizeBytes',0);
entries = repmat(template,numel(selected),1);
for k = 1:numel(selected)
    entries(k).name = [char(string(selected(k).sceneName)) '.zip'];
    entries(k).url = char(string(selected(k).url));
    if isfield(selected,'sizeBytes')
        entries(k).sizeBytes = double(selected(k).sizeBytes);
    end
end
end

function cleanupTransferDirectory(pathValue)
if isempty(pathValue) || ~isfolder(pathValue), return; end
try, rmdir(pathValue,'s'); catch, end
end

function openFolder(pathValue)
if ispc, winopen(pathValue);
else, desktop = java.awt.Desktop.getDesktop(); desktop.open(java.io.File(pathValue)); end
end

function polygon = bboxPolygon(cfg)
polygon = [cfg.lon_min cfg.lat_min; cfg.lon_max cfg.lat_min; ...
    cfg.lon_max cfg.lat_max; cfg.lon_min cfg.lat_max; cfg.lon_min cfg.lat_min];
end

function value = numericMatrix(raw)
if iscell(raw)
    try
        if isvector(raw) && all(cellfun(@(row) isnumeric(row) && numel(row) == 2, raw))
            value = vertcat(raw{:});
        else
            value = cell2mat(raw);
        end
    catch
        value = zeros(0, 2);
    end
else
    value = double(raw);
end
value = squeeze(value);
if size(value, 2) ~= 2 && size(value, 1) == 2, value = value.'; end
end

function changed = bboxChanged(left, right)
fields = {'lon_min','lon_max','lat_min','lat_max'};
changed = any(cellfun(@(name) ~isequaln(left.(name), right.(name)), fields));
end

function closeProgress(dialog)
try, if ~isempty(dialog) && isvalid(dialog), close(dialog); end, catch, end
end

function filters = normalizeFilters(raw, filters)
multiFields = {'processingLevel','beamMode','polarization','flightDirection','subtype'};
for k = 1:numel(multiFields)
    name = multiFields{k};
    if isfield(raw,name), filters.(name) = asCellText(raw.(name)); end
end
textFields = {'dataset','startDate','endDate','pathStart','pathEnd', ...
    'frameStart','frameEnd','groupID','samplingRate','samplingUnit'};
for k = 1:numel(textFields)
    name = textFields{k};
    if isfield(raw,name), filters.(name) = char(string(raw.(name))); end
end
filters.dataset = 'SENTINEL-1';
end

function request = buildSearchRequest(filters, polygon)
request = struct();
request.repository = 'ASF'; request.dataset = 'SENTINEL-1';
request.processingLevel = filters.processingLevel;
request.beamMode = filters.beamMode;
request.polarization = filters.polarization;
request.flightDirection = filters.flightDirection;
request.subtype = filters.subtype;
request.startDate = optionalText(filters.startDate);
request.endDate = optionalText(filters.endDate);
request.pathStart = optionalNumber(filters.pathStart);
request.pathEnd = optionalNumber(filters.pathEnd);
request.frameStart = optionalNumber(filters.frameStart);
request.frameEnd = optionalNumber(filters.frameEnd);
request.groupID = optionalNumber(filters.groupID);
request.sampling = struct('rate',optionalNumber(filters.samplingRate),'unit',[]);
if ~isempty(request.sampling.rate), request.sampling.unit = char(string(filters.samplingUnit)); end
corners = polygon;
if size(corners,1)>1 && isequal(corners(1,:),corners(end,:)), corners(end,:) = []; end
request.aoi = struct('type','polygon','corners',corners);
end

function value = optionalText(raw)
value = strtrim(char(string(raw)));
if isempty(value), value = []; end
end

function value = optionalNumber(raw)
if isnumeric(raw) && isscalar(raw), value = double(raw); else, value = str2double(strtrim(char(string(raw)))); end
if isempty(value) || ~isscalar(value) || ~isfinite(value) || value==0, value = []; end
end

function values = asCellText(raw)
if isempty(raw), values = {}; return; end
if iscell(raw), values = cellfun(@(item) char(string(item)),raw,'UniformOutput',false);
else, values = cellstr(string(raw)); end
values = values(~cellfun(@isempty,values));
end

function writeJson(pathValue, value)
folder = fileparts(pathValue); if ~isfolder(folder), mkdir(folder); end
fid = fopen(pathValue,'w');
if fid==-1, error('PHASE:CannotWriteJSON','Cannot write %s.',pathValue); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end

function backupAsfFiles(rootDir)
folder = fullfile(rootDir,'downloadasf');
pairs = {'download_data.json','download_data_all.json'; ...
    'search_summary.json','search_summary_all.json'};
for k = 1:size(pairs,1)
    source = fullfile(folder,pairs{k,1}); target = fullfile(folder,pairs{k,2});
    if isfile(source), copyfile(source,target,'f'); end
end
end

function names = payloadNames(payload)
if ~isstruct(payload) || ~isfield(payload,'sceneNames'), names = strings(0,1); return; end
names = string(payload.sceneNames); names = names(strlength(names)>0);
end

function validateCompatibleSelection(products)
if isempty(products), error('PHASE:NoASFSelection','Select at least one ASF result.'); end
paths = unique([products.pathNumber]); frames = unique([products.frameNumber]);
directions = unique(string({products.direction}));
if numel(paths)>1 || numel(frames)>1 || numel(directions)>1
    error('PHASE:IncompatibleASFSelection', ...
        'Selected products must share one path, frame and direction. Paths: %s; frames: %s; directions: %s.', ...
        mat2str(paths),mat2str(frames),strjoin(directions,', '));
end
end

function count = countExisting(folder, products)
count = 0;
for k = 1:numel(products)
    if isfile(fullfile(folder,[products(k).sceneName '.zip'])), count = count+1; end
end
end

function saveAsfSelection(rootDir, products)
folder = fullfile(rootDir,'downloadasf');
downloadPath = fullfile(folder,'download_data_all.json');
summaryPath = fullfile(folder,'search_summary_all.json');
if ~isfile(downloadPath) || ~isfile(summaryPath)
    error('PHASE:MissingASFSearchFiles','Run an ASF search before downloading.');
end
names = string({products.sceneName});
download = jsondecode(fileread(downloadPath));
information = download.information;
selectedInfo = information(ismember(string({information.sceneName}),names));
output = struct(); output.information = selectedInfo;
writeJson(fullfile(folder,'download_data.json'),output);

summary = jsondecode(fileread(summaryPath));
allProducts = summary.products;
selectedProducts = allProducts(ismember(string({allProducts.sceneName}),names));
summary.products = selectedProducts;
summary.product_count = numel(selectedProducts);
summary.total_size_gb = round(sum([selectedProducts.size]),2);
summary.total_size_bytes = sum([selectedProducts.size_bytes]);
writeJson(fullfile(folder,'search_summary.json'),summary);

requestPath = fullfile(folder,'search_request.json');
if isfile(requestPath)
    request = jsondecode(fileread(requestPath));
    request.pathStart = products(1).pathNumber; request.pathEnd = products(1).pathNumber;
    request.frameStart = products(1).frameNumber; request.frameEnd = products(1).frameNumber;
    request.flightDirection = {products(1).direction};
    writeJson(fullfile(folder,'last_download_request.json'),request);
end
end

function value = quoteCommandArgument(raw)
value = char(string(raw)); value = strrep(value,'"','""'); value = ['"' value '"'];
end

function deleteIfExists(pathValue)
if isfile(pathValue), delete(pathValue); end
end

function results = normalizeUpdateResults(raw)
template = struct('sceneName','','date','','startTime','','platform','', ...
    'polarization','','pathNumber',0,'frameNumber',0, ...
    'url','','sizeBytes',0,'selected',true);
if isempty(raw), results = repmat(template,0,1); return; end
results = repmat(template,numel(raw),1);
for k = 1:numel(raw)
    fields = {'sceneName','date','startTime','platform','polarization','url'};
    for j = 1:numel(fields)
        name = fields{j};
        if isfield(raw(k),name), results(k).(name) = char(string(raw(k).(name))); end
    end
    if isfield(raw(k),'pathNumber'), results(k).pathNumber = double(raw(k).pathNumber); end
    if isfield(raw(k),'frameNumber'), results(k).frameNumber = double(raw(k).frameNumber); end
    if isfield(raw(k),'sizeBytes'), results(k).sizeBytes = double(raw(k).sizeBytes); end
    results(k).selected = true;
end
end

function writeLines(pathValue, lines)
fid = fopen(pathValue,'w');
if fid==-1, error('PHASE:CannotWriteFile','Cannot write %s.',pathValue); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
for k = 1:numel(lines), fprintf(fid,'%s\n',char(string(lines{k}))); end
end

function cleanupFiles(paths)
for k = 1:numel(paths), deleteIfExists(paths{k}); end
end

function keys = successfulTileKeys(result)
keys = {};
if ~isstruct(result) || ~isfield(result,'successful') || isempty(result.successful)
    return;
end
successful = result.successful;
if isstruct(successful) && isfield(successful,'key')
    keys = cellstr(string({successful.key}));
end
end
