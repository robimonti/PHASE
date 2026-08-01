classdef App < handle
    %APP Liquid Glass controller around the editable PHASE Model engine.

    properties (SetAccess = private)
        UIFigure
        HTML
        Engine
        RootDir
        Config
        SavedConfig
        Logs = {}
        Status = 'idle'
        StatusDetail = 'Initialising PHASE Model'
        IsDirty = false
        IsRunning = false
        StopRequested = false
        RunStartedAt = []
        Progress = struct('percentage',0,'phase','Ready', ...
            'elapsedSeconds',0)
        LiveLogFile = ''
        LiveLogUrl = ''
        DiaryActive = false
        MapBase
        MapPolygon
        MapPsExtent = zeros(0,2)
        MapAoiFootprints = struct([])
        MapTileBusy = false
    end

    methods
        function obj = App(rootDir)
            if nargin < 1 || isempty(rootDir)
                rootDir = phase_model_beta.projectRoot();
            end
            obj.RootDir = char(java.io.File(rootDir).getCanonicalPath());
            [obj.Config,info] = phase_model_beta.loadConfig(obj.RootDir);
            obj.SavedConfig = obj.Config;
            obj.IsDirty = ~info.exists;
            obj.MapBase = phase_model_beta.mapBase();
            obj.MapPolygon = configPolygon(obj.Config);
            obj.applyDefaultPsBounds(false);
            obj.refreshMapAoi(false);
            obj.configureLiveLog();
            obj.ensureAssets();

            obj.UIFigure = uifigure( ...
                'Name','PHASE · Geospatial Model', ...
                'Color',[1 1 1], ...
                'Position',centeredPosition(1500,920));
            obj.UIFigure.UserData = obj;
            obj.UIFigure.CloseRequestFcn = @(~,~) delete(obj);
            grid = uigridlayout(obj.UIFigure,[1 1]);
            grid.Padding = [0 0 0 0];
            uiPath = fullfile(obj.RootDir,'phase_model_beta_ui','index.html');
            obj.HTML = uihtml(grid,'HTMLSource',uiPath);
            obj.HTML.HTMLEventReceivedFcn = @(~,event) obj.onHtmlEvent(event);
            drawnow;
            obj.sendState();

            try
                obj.appendLog('Initialising the standalone PHASE Model engine…');
                obj.Engine = phase_model_beta.LegacyEngine();
                obj.Engine.ExternalLogCallback = @(message) obj.appendLog(message);
                obj.Engine.ExternalProgressCallback = @(progress) obj.onProgress(progress);
                obj.Engine.UIFigure.CloseRequestFcn = @(~,~) obj.hideEngine();
                obj.Engine.UIFigure.Visible = 'off';
                phase_model_beta.applyConfigToEngine(obj.Engine,obj.Config);
                if info.exists && ~obj.IsDirty
                    obj.Status = 'saved';
                    obj.StatusDetail = 'Configuration loaded; model ready';
                    obj.appendLog(['Loaded configuration: ' info.path]);
                elseif info.exists
                    obj.Status = 'idle';
                    obj.StatusDetail = 'Default PS extent detected; save before starting';
                    obj.appendLog(['Loaded configuration and initialised its missing AOI: ' ...
                        info.path]);
                else
                    obj.Status = 'idle';
                    obj.StatusDetail = 'Review and save the initial configuration';
                    obj.appendLog('No input_model.mat exists yet; initial defaults loaded.');
                end
                obj.appendLog('PHASE Model engine ready.');
            catch ME
                obj.Status = 'error';
                obj.StatusDetail = ME.message;
                obj.appendLog(['Engine initialisation failed [' ME.identifier ']: ' ME.message]);
                obj.showError('PHASE Model',ME.message);
            end
            obj.sendState();
        end

        function appendLog(obj,message)
            message = char(string(message));
            if isempty(strtrim(message)), return; end
            entry = struct( ...
                'time',char(datetime('now','Format','HH:mm:ss')), ...
                'message',message);
            obj.Logs{end+1} = entry;
            if numel(obj.Logs) > 800
                obj.Logs = obj.Logs(end-799:end);
            end
            try
                if ~isempty(obj.HTML) && isvalid(obj.HTML)
                    sendEventToHTMLSource(obj.HTML,'PhaseLog',entry);
                    drawnow limitrate
                end
            catch
            end
        end

        function delete(obj)
            obj.endDiary();
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
        function onHtmlEvent(obj,event)
            try
                name = lower(char(string(event.HTMLEventName)));
                payload = event.HTMLEventData;
                switch name
                    case 'ready'
                        obj.sendState();
                    case 'changed'
                        obj.updateFromPayload(payload);
                    case 'load'
                        obj.loadFromDisk();
                    case 'save'
                        obj.saveFromPayload(payload);
                    case 'start'
                        obj.startFromPayload(payload);
                    case 'stop'
                        obj.requestHardStop();
                    case 'browse'
                        obj.browseFromPayload(payload);
                    case 'estimatepsbounds'
                        obj.estimatePsBoundsFromPayload(payload);
                    case 'mapaoichanged'
                        obj.mapAoiChanged(payload);
                    case 'maptilesrequested'
                        obj.cacheMapTiles(payload);
                    case 'openroot'
                        openFolder(obj.RootDir);
                    case 'clearlog'
                        obj.Logs = {};
                        obj.sendState();
                    otherwise
                        obj.appendLog(['Unknown interface event: ' name]);
                end
            catch ME
                if ~obj.IsRunning
                    obj.Status = 'error';
                    obj.StatusDetail = ME.message;
                end
                obj.appendLog(['Interface action failed [' ME.identifier ']: ' ME.message]);
                obj.showError('PHASE Model',ME.message);
                obj.sendState();
            end
        end

        function updateFromPayload(obj,payload)
            candidate = obj.payloadConfig(payload);
            inputChanged = ~strcmp(char(string(candidate.filepathIN)), ...
                char(string(obj.Config.filepathIN)));
            obj.Config = candidate;
            if inputChanged
                obj.applyDefaultPsBounds(true);
            end
            obj.refreshMapAoi(false);
            obj.IsDirty = ~phase_model_beta.configsEqual(obj.Config,obj.SavedConfig);
            if ~obj.IsRunning
                obj.Status = ternary(obj.IsDirty,'idle','saved');
                obj.StatusDetail = ternary(obj.IsDirty, ...
                    'Unsaved changes','Configuration saved and ready');
            end
            obj.sendState();
        end

        function saveFromPayload(obj,payload)
            if obj.IsRunning, return; end
            candidate = obj.payloadConfig(payload);
            [errors,warnings] = phase_model_beta.validateConfig(candidate,false);
            if ~isempty(errors)
                error('PHASE_Model_beta:invalidConfiguration','%s',strjoin(errors,newline));
            end
            pathValue = phase_model_beta.saveConfig(obj.RootDir,candidate);
            obj.Config = candidate;
            obj.SavedConfig = candidate;
            obj.refreshMapAoi(true);
            obj.IsDirty = false;
            phase_model_beta.applyConfigToEngine(obj.Engine,candidate);
            obj.Engine.UIFigure.Visible = 'off';
            obj.Status = 'saved';
            obj.StatusDetail = 'Configuration saved and ready';
            obj.appendLog(['Saved configuration: ' pathValue]);
            for k = 1:numel(warnings), obj.appendLog(['NOTICE: ' warnings{k}]); end
            obj.sendState();
        end

        function loadFromDisk(obj)
            if obj.IsRunning, return; end
            [candidate,info] = phase_model_beta.loadConfig(obj.RootDir);
            if ~info.exists
                error('PHASE_Model_beta:configMissing', ...
                    'No input_model.mat exists yet. Review the settings and press Save.');
            end
            obj.Config = candidate;
            obj.SavedConfig = candidate;
            obj.MapPolygon = configPolygon(candidate);
            obj.IsDirty = false;
            obj.applyDefaultPsBounds(false);
            obj.refreshMapAoi(true);
            phase_model_beta.applyConfigToEngine(obj.Engine,obj.Config);
            obj.Engine.UIFigure.Visible = 'off';
            obj.Status = ternary(obj.IsDirty,'idle','saved');
            obj.StatusDetail = ternary(obj.IsDirty, ...
                'Default PS extent detected; save before starting', ...
                'Configuration loaded');
            obj.appendLog(['Loaded configuration: ' info.path]);
            obj.sendState();
        end

        function startFromPayload(obj,payload)
            if obj.IsRunning, return; end
            candidate = obj.payloadConfig(payload);
            [errors,warnings] = phase_model_beta.validateConfig(candidate,true);
            if ~isempty(errors)
                error('PHASE_Model_beta:invalidConfiguration','%s',strjoin(errors,newline));
            end
            if ~phase_model_beta.configsEqual(candidate,obj.SavedConfig)
                error('PHASE_Model_beta:unsavedConfiguration', ...
                    'Save the visible configuration before starting.');
            end
            obj.Config = candidate;
            obj.IsDirty = false;
            obj.IsRunning = true;
            obj.StopRequested = false;
            obj.Engine.StopRequested = false;
            obj.Status = 'running';
            obj.StatusDetail = 'Preparing environment';
            obj.RunStartedAt = datetime('now');
            obj.Progress = struct('percentage',0,'phase','Preparing environment', ...
                'elapsedSeconds',0);
            obj.beginDiary();
            obj.appendLog('Starting PHASE geospatial modelling.');
            for k = 1:numel(warnings), obj.appendLog(['NOTICE: ' warnings{k}]); end
            obj.sendState();
            try
                phase_model_beta.applyConfigToEngine(obj.Engine,candidate);
                obj.Engine.UIFigure.Visible = 'off';
                obj.Engine.StartButtonPushed([]);
                obj.Engine.UIFigure.Visible = 'off';
                obj.IsRunning = false;
                obj.Status = 'success';
                if isempty(obj.Engine.outputDir)
                    obj.StatusDetail = 'Processing completed';
                else
                    obj.StatusDetail = ['Completed: ' obj.Engine.outputDir];
                end
                obj.Progress.percentage = 100;
                obj.Progress.phase = 'Processing complete';
                obj.Progress.elapsedSeconds = elapsedSeconds(obj.RunStartedAt);
                obj.appendLog(obj.StatusDetail);
            catch ME
                obj.IsRunning = false;
                if strcmp(ME.identifier,'PHASE_Model_beta:hardStopped')
                    obj.Status = 'stopped';
                    obj.StatusDetail = 'Processing stopped by the user';
                    obj.Progress.phase = 'Stopped';
                    obj.appendLog('Hard stop completed at a safe numerical checkpoint.');
                else
                    obj.Status = 'error';
                    obj.StatusDetail = ME.message;
                    obj.Progress.phase = 'Processing failed';
                    obj.appendLog(['Processing failed [' ME.identifier ']: ' ME.message]);
                    obj.appendLog(getReport(ME,'extended','hyperlinks','off'));
                    obj.showError('PHASE Model processing failed',ME.message);
                end
                obj.Progress.elapsedSeconds = elapsedSeconds(obj.RunStartedAt);
            end
            obj.Engine.StopRequested = false;
            obj.StopRequested = false;
            obj.endDiary();
            obj.sendState();
        end

        function browseFromPayload(obj,payload)
            if obj.IsRunning, return; end
            if isstruct(payload) && isfield(payload,'config')
                obj.updateFromPayload(payload.config);
            end
            fieldName = char(string(payload.field));
            selected = '';
            switch fieldName
                case 'filepathIN'
                    [file,path] = uigetfile({'*.xlsx;*.csv','Time-series files (*.xlsx, *.csv)'}, ...
                        'Select displacement time series');
                    if ~isequal(file,0), selected = fullfile(path,file); end
                case 'pythonPath'
                    [file,path] = uigetfile({'*','Python executable'}, ...
                        'Select Python executable');
                    if ~isequal(file,0), selected = fullfile(path,file); end
                case 'filepathAOI'
                    [file,path] = uigetfile({'*.shp','Shapefile (*.shp)'}, ...
                        'Select AOI shapefile');
                    if ~isequal(file,0), selected = fullfile(path,file); end
                case 'filepath_EXTR'
                    [file,path] = uigetfile({'*.txt;*.csv','Query-point files (*.txt, *.csv)'}, ...
                        'Select query-points file');
                    if ~isequal(file,0), selected = fullfile(path,file); end
                case 'coherence_dir'
                    value = uigetdir(obj.RootDir,'Select coherence folder');
                    if ~isequal(value,0), selected = value; end
                otherwise
                    error('PHASE_Model_beta:invalidBrowseField', ...
                        'Unsupported browse field: %s',fieldName);
            end
            if isempty(selected), return; end
            obj.Config.(fieldName) = selected;
            if strcmp(fieldName,'filepathIN')
                obj.applyDefaultPsBounds(true);
            elseif strcmp(fieldName,'filepathAOI')
                obj.Config.flag_AOIbb = false;
            end
            obj.refreshMapAoi(true);
            obj.IsDirty = ~phase_model_beta.configsEqual(obj.Config,obj.SavedConfig);
            obj.Status = 'idle';
            obj.StatusDetail = 'Unsaved changes';
            obj.sendState();
        end

        function estimatePsBoundsFromPayload(obj,payload)
            if obj.IsRunning, return; end
            candidate = obj.payloadConfig(payload);
            bounds = phase_model_beta.estimatePsBoundingBox(candidate.filepathIN);
            candidate.flag_AOIbb = true;
            candidate.lonMinAOI = bounds(1);
            candidate.lonMaxAOI = bounds(2);
            candidate.latMinAOI = bounds(3);
            candidate.latMaxAOI = bounds(4);
            candidate.aoi_polygon_lonlat = bboxPolygon(bounds);
            obj.Config = candidate;
            obj.MapPolygon = candidate.aoi_polygon_lonlat;
            obj.MapAoiFootprints = struct([]);
            obj.IsDirty = ~phase_model_beta.configsEqual(candidate,obj.SavedConfig);
            obj.Status = 'idle';
            obj.StatusDetail = 'Full PS extent estimated; save before starting';
            obj.appendLog(sprintf( ...
                'AOI set to full PS extent: lon %.6f to %.6f, lat %.6f to %.6f.', ...
                bounds(1),bounds(2),bounds(3),bounds(4)));
            obj.sendState();
        end

        function applyDefaultPsBounds(obj,force)
            if nargin < 2, force = false; end
            inputPath = char(string(obj.Config.filepathIN));
            obj.MapPsExtent = zeros(0,2);
            if isempty(strtrim(inputPath)) || ~isfile(inputPath)
                obj.MapPolygon = configPolygon(obj.Config);
                return
            end
            try
                bounds = phase_model_beta.estimatePsBoundingBox(inputPath);
                obj.MapPsExtent = bboxPolygon(bounds);
                existingPolygon = configPolygon(obj.Config);
                hasShapefile = ~obj.Config.flag_AOIbb && ...
                    ~isempty(strtrim(char(string(obj.Config.filepathAOI))));
                if ~isempty(existingPolygon) || hasShapefile
                    obj.MapPolygon = existingPolygon;
                    return
                end
                obj.Config.flag_AOIbb = true;
                obj.Config.lonMinAOI = bounds(1);
                obj.Config.lonMaxAOI = bounds(2);
                obj.Config.latMinAOI = bounds(3);
                obj.Config.latMaxAOI = bounds(4);
                obj.Config.aoi_polygon_lonlat = bboxPolygon(bounds);
                obj.MapPolygon = obj.Config.aoi_polygon_lonlat;
                obj.IsDirty = ~phase_model_beta.configsEqual( ...
                    obj.Config,obj.SavedConfig);
                obj.appendLog(sprintf( ...
                    ['Default AOI fitted to the PS extent: lon %.6f to %.6f, ' ...
                     'lat %.6f to %.6f.'], ...
                    bounds(1),bounds(2),bounds(3),bounds(4)));
            catch ME
                obj.appendLog(['Could not initialise the AOI from the PS file [' ...
                    ME.identifier ']: ' ME.message]);
            end
        end

        function refreshMapAoi(obj,logFailures)
            if nargin < 2, logFailures = false; end
            obj.MapAoiFootprints = struct([]);
            if obj.Config.flag_AOIbb
                obj.MapPolygon = configPolygon(obj.Config);
                return
            end
            obj.MapPolygon = zeros(0,2);
            shapefilePath = char(string(obj.Config.filepathAOI));
            if isempty(strtrim(shapefilePath)) || ~isfile(shapefilePath), return; end
            try
                [~,segments,info] = phase_model_beta.readAoiShapefile( ...
                    shapefilePath,obj.Config.filepathIN);
                [~,shapefileName,extension] = fileparts(shapefilePath);
                features = repmat(struct('id','','name','','source','', ...
                    'selected',true,'coordinates',zeros(0,2)),1,numel(segments));
                for partIndex = 1:numel(segments)
                    features(partIndex).id = sprintf('aoi-shapefile-%d',partIndex);
                    features(partIndex).name = sprintf('Selected AOI · part %d',partIndex);
                    features(partIndex).source = [shapefileName extension];
                    features(partIndex).selected = true;
                    features(partIndex).coordinates = segments{partIndex};
                end
                obj.MapAoiFootprints = features;
                if logFailures
                    obj.appendLog(sprintf( ...
                        'AOI shapefile displayed on map: %s (%d polygon part(s), %s coordinates).', ...
                        shapefilePath,info.partCount,info.coordinateType));
                end
            catch ME
                if logFailures
                    obj.appendLog(['Could not display the selected AOI shapefile [' ...
                        ME.identifier ']: ' ME.message]);
                end
            end
        end

        function requestHardStop(obj)
            if ~obj.IsRunning || obj.StopRequested, return; end
            obj.StopRequested = true;
            if ~isempty(obj.Engine) && isvalid(obj.Engine)
                obj.Engine.StopRequested = true;
            end
            obj.StatusDetail = 'Hard stop requested; stopping at the next safe checkpoint';
            obj.Progress.phase = 'Stopping…';
            obj.appendLog('HARD STOP requested by the user.');
            obj.sendStatus();
            drawnow;
        end

        function mapAoiChanged(obj,payload)
            if obj.IsRunning, return; end
            if ~isstruct(payload) || ~isfield(payload,'polygon') || ~isfield(payload,'bbox')
                error('PHASE_Model_beta:invalidMapAOI', ...
                    'The map did not provide polygon and bounding-box data.');
            end
            polygon = numericMatrix(payload.polygon);
            bbox = payload.bbox;
            required = {'minLon','maxLon','minLat','maxLat'};
            if size(polygon,1)<3 || size(polygon,2)~=2 || ...
                    ~all(cellfun(@(name) isfield(bbox,name),required))
                error('PHASE_Model_beta:invalidMapAOI', ...
                    'Draw a polygon with at least three valid vertices.');
            end
            values = cellfun(@(name) double(bbox.(name)),required);
            if any(~isfinite(values)) || values(1)>=values(2) || values(3)>=values(4)
                error('PHASE_Model_beta:invalidMapBounds','The polygon bounding box is invalid.');
            end
            if ~isequal(polygon(1,:),polygon(end,:)), polygon(end+1,:)=polygon(1,:); end
            obj.MapPolygon = polygon;
            obj.MapAoiFootprints = struct([]);
            obj.Config.flag_AOIbb = true;
            obj.Config.aoi_polygon_lonlat = polygon;
            obj.Config.lonMinAOI = values(1); obj.Config.lonMaxAOI = values(2);
            obj.Config.latMinAOI = values(3); obj.Config.latMaxAOI = values(4);
            obj.IsDirty = ~phase_model_beta.configsEqual(obj.Config,obj.SavedConfig);
            obj.Status = 'idle'; obj.StatusDetail = 'AOI polygon updated; save before starting';
            obj.appendLog(sprintf( ...
                'AOI polygon updated (%d vertices): lon %.6f..%.6f, lat %.6f..%.6f.', ...
                size(polygon,1)-1,values(1),values(2),values(3),values(4)));
            obj.sendState();
        end

        function cacheMapTiles(obj,payload)
            if obj.MapTileBusy || ~isstruct(payload) || ~isfield(payload,'requests'), return; end
            requests = payload.requests;
            if isempty(requests), return; end
            if numel(requests)>120
                obj.appendLog('Map tile request ignored because it exceeded the safe batch limit.');
                return
            end
            obj.MapTileBusy = true;
            requestPath = [tempname '.json']; outputPath = [tempname '.json'];
            cleanup = onCleanup(@() cleanupFiles({requestPath,outputPath})); %#ok<NASGU>
            try
                writeJson(requestPath,struct('requests',{requests}));
                cacheRoot = fullfile(obj.RootDir,'phase_model_beta_ui','map_tiles');
                if ~isfolder(cacheRoot), mkdir(cacheRoot); end
                script = fullfile(obj.RootDir,'pythonScripts','cache_phase_map_tiles.py');
                parts = {phase_model_beta.resolvePythonPath(obj.Config.pythonPath),script, ...
                    '--cache-root',cacheRoot,'--requests',requestPath,'--output',outputPath};
                quoted = cellfun(@quoteCommandArgument,parts,'UniformOutput',false);
                [status,output] = phase_model_beta.runCommandHidden( ...
                    obj.Engine,strjoin(quoted,' '),'Satellite map tiles');
                if status~=0 || ~isfile(outputPath)
                    error('PHASE:MapTileCacheFailed','%s',strtrim(output));
                end
                result = jsondecode(fileread(outputPath));
                keys = successfulTileKeys(result); failedCount = 0;
                if isfield(result,'failed'), failedCount = numel(result.failed); end
                sendEventToHTMLSource(obj.HTML,'MapTilesReady', ...
                    struct('keys',{keys},'failedCount',failedCount));
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

        function config = payloadConfig(~,payload)
            if isstruct(payload) && isfield(payload,'config')
                payload = payload.config;
            end
            if ~isstruct(payload)
                error('PHASE_Model_beta:missingConfig', ...
                    'The interface did not send a configuration.');
            end
            config = phase_model_beta.uiToConfig(payload);
        end

        function onProgress(obj,progress)
            obj.Progress.percentage = max(0,min(100,double(progress.percentage)));
            obj.Progress.phase = char(string(progress.phase));
            obj.Progress.elapsedSeconds = elapsedSeconds(obj.RunStartedAt);
            obj.StatusDetail = obj.Progress.phase;
            obj.sendStatus();
        end

        function sendState(obj)
            if isempty(obj.HTML) || ~isvalid(obj.HTML), return; end
            if isempty(obj.MapPsExtent)
                mapFootprints = struct([]);
            else
                mapFootprints = struct( ...
                    'id','ps-extent','name','PS extent', ...
                    'source','Selected displacement file','selected',false, ...
                    'coordinates',obj.MapPsExtent);
            end
            if ~isempty(obj.MapAoiFootprints)
                if isempty(mapFootprints)
                    mapFootprints = obj.MapAoiFootprints;
                else
                    mapFootprints = [mapFootprints obj.MapAoiFootprints];
                end
            end
            mapState = struct('coastlines',{obj.MapBase}, ...
                'footprints',{mapFootprints},'polygon',obj.MapPolygon);
            state = struct( ...
                'kind','state', ...
                'version','6.0.0', ...
                'schema',phase_model_beta.schema(), ...
                'config',phase_model_beta.configToUi(obj.Config), ...
                'rootDir',obj.RootDir, ...
                'status',obj.Status, ...
                'statusDetail',obj.StatusDetail, ...
                'dirty',obj.IsDirty, ...
                'running',obj.IsRunning, ...
                'progress',obj.Progress, ...
                'map',mapState, ...
                'logs',{obj.Logs}, ...
                'liveLogUrl',obj.LiveLogUrl);
            try, obj.HTML.Data = state; drawnow limitrate; catch, end
        end

        function sendStatus(obj)
            if isempty(obj.HTML) || ~isvalid(obj.HTML), return; end
            state = struct( ...
                'status',obj.Status, ...
                'statusDetail',obj.StatusDetail, ...
                'dirty',obj.IsDirty, ...
                'running',obj.IsRunning, ...
                'progress',obj.Progress);
            try
                sendEventToHTMLSource(obj.HTML,'PhaseStatus',state);
                drawnow limitrate
            catch
            end
        end

        function configureLiveLog(obj)
            runtimeDir = fullfile(obj.RootDir,'phase_model_beta_ui','runtime_logs');
            if ~isfolder(runtimeDir), mkdir(runtimeDir); end
            fileName = ['model_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')) '.log'];
            obj.LiveLogFile = fullfile(runtimeDir,fileName);
            obj.LiveLogUrl = ['runtime_logs/' fileName];
        end

        function beginDiary(obj)
            obj.endDiary();
            try
                if isfile(obj.LiveLogFile), delete(obj.LiveLogFile); end
                diary(obj.LiveLogFile);
                diary on
                obj.DiaryActive = true;
            catch ME
                obj.appendLog(['Live diary unavailable: ' ME.message]);
            end
        end

        function endDiary(obj)
            if ~obj.DiaryActive, return; end
            try, diary off; catch, end
            obj.DiaryActive = false;
        end

        function ensureAssets(obj)
            assetsDir = fullfile(obj.RootDir,'phase_model_beta_ui','assets');
            if ~isfolder(assetsDir), mkdir(assetsDir); end
            assets = {
                fullfile(obj.RootDir,'PHASE_logo.png'), fullfile(assetsDir,'PHASE_logo.png')
                fullfile(obj.RootDir,'PHASE_mod2.png'), fullfile(assetsDir,'PHASE_mod2.png')
            };
            for k = 1:size(assets,1)
                if isfile(assets{k,1}) && ~isfile(assets{k,2})
                    copyfile(assets{k,1},assets{k,2});
                end
            end

            % uihtml on Windows does not reliably allow a page to load a
            % script through "../" from a sibling directory. Keep the shared
            % map implementation authoritative, but stage a local copy beside
            % the Model page before the HTML component is constructed.
            mapSource = fullfile(obj.RootDir,'PHASE_Preprocessing', ...
                'phase_preprocessing_beta_ui','map.js');
            mapTarget = fullfile(obj.RootDir,'phase_model_beta_ui','map.js');
            if ~isfile(mapSource)
                error('PHASE_Model_beta:mapRuntimeMissing', ...
                    'The shared PHASE map runtime is missing: %s', mapSource);
            end
            [copied,message] = copyfile(mapSource,mapTarget,'f');
            if ~copied
                error('PHASE_Model_beta:mapRuntimeCopyFailed', ...
                    'Could not prepare the local Model map runtime: %s', message);
            end
        end

        function hideEngine(obj)
            try, obj.Engine.UIFigure.Visible = 'off'; catch, end
        end

        function showError(obj,titleText,message)
            try
                uialert(obj.UIFigure,message,titleText,'Icon','error');
            catch
            end
        end
    end
end

function position = centeredPosition(width,height)
screen = get(groot,'ScreenSize');
width = min(width,max(1000,screen(3)-80));
height = min(height,max(720,screen(4)-100));
position = [max(1,(screen(3)-width)/2) max(1,(screen(4)-height)/2) width height];
end

function value = ternary(condition,whenTrue,whenFalse)
if condition, value = whenTrue; else, value = whenFalse; end
end

function value = elapsedSeconds(startedAt)
if isempty(startedAt), value = 0; else, value = max(0,seconds(datetime('now')-startedAt)); end
end

function openFolder(pathValue)
if ispc
    winopen(pathValue);
elseif ismac
    system(['open ' quoteShell(pathValue)]);
else
    system(['xdg-open ' quoteShell(pathValue) ' >/dev/null 2>&1 &']);
end
end

function value = quoteShell(raw)
value = ['''' strrep(char(string(raw)),'''','''"''"''') ''''];
end

function polygon = configPolygon(config)
polygon = config.aoi_polygon_lonlat;
if ~config.flag_AOIbb
    polygon = zeros(0,2);
elseif isempty(polygon)
    bounds = [config.lonMinAOI config.lonMaxAOI config.latMinAOI config.latMaxAOI];
    if all(isfinite(bounds)) && bounds(1)<bounds(2) && bounds(3)<bounds(4)
        polygon = bboxPolygon(bounds);
    else
        polygon = zeros(0,2);
    end
end
end

function polygon = bboxPolygon(bounds)
polygon = [bounds(1) bounds(3); bounds(2) bounds(3); ...
    bounds(2) bounds(4); bounds(1) bounds(4); bounds(1) bounds(3)];
end

function value = numericMatrix(raw)
if iscell(raw)
    try, value = cell2mat(raw); catch, value = zeros(0,2); end
else
    value = double(raw);
end
value = squeeze(value);
if size(value,2)~=2 && size(value,1)==2, value=value.'; end
end

function writeJson(pathValue,value)
fid = fopen(pathValue,'w');
if fid==-1, error('PHASE:CannotWriteJSON','Cannot write %s.',pathValue); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value));
end

function value = quoteCommandArgument(raw)
value = char(string(raw)); value = strrep(value,'"','""'); value = ['"' value '"'];
end

function cleanupFiles(paths)
for k=1:numel(paths)
    if isfile(paths{k}), try, delete(paths{k}); catch, end, end
end
end

function keys = successfulTileKeys(result)
keys = {};
if isstruct(result) && isfield(result,'successful') && ...
        isstruct(result.successful) && isfield(result.successful,'key')
    keys = cellstr(string({result.successful.key}));
end
end
