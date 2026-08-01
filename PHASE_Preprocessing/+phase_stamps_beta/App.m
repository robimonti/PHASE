classdef App < handle
    %APP Controller for the HTML/CSS/JavaScript PHASE StaMPS beta.

    properties (SetAccess = private)
        UIFigure
        HTML
        WorkDir
        LauncherDir
        Config
        SavedConfig
        AutoDetectedFields = {}
        Logs = {}
        Status = 'idle'
        StatusDetail = 'Ready'
        IsDirty = false
        IsRunning = false
        TSPickerOverlay = []
        TSPickerContainer = []
        LiveLogFile = ''
        LiveLogUrl = ''
        DiaryActive = false
    end

    methods
        function obj = App(workDir, launcherDir)
            obj.WorkDir = char(java.io.File(workDir).getCanonicalPath());
            obj.LauncherDir = launcherDir;

            [loaded, info] = phase_stamps_beta.loadConfig(obj.WorkDir);
            obj.SavedConfig = loaded;
            [detected, fields, messages] = ...
                phase_stamps_beta.autoDetectConfig(loaded, obj.WorkDir);
            obj.Config = detected;
            obj.AutoDetectedFields = fields;
            obj.IsDirty = ~phase_stamps_beta.configsEqual( ...
                obj.Config, obj.SavedConfig) || ~info.exists;

            obj.UIFigure = uifigure('Name', 'PHASE · StaMPS', ...
                'Color', [0.025 0.039 0.075], ...
                'Position', centeredPosition(1440, 900));
            % Keep the controller alive when the launcher is called without
            % an output argument. The figure releases this reference on close.
            obj.UIFigure.UserData = obj;
            obj.UIFigure.CloseRequestFcn = @(~,~) delete(obj);
            grid = uigridlayout(obj.UIFigure, [1 1]);
            grid.Padding = [0 0 0 0];
            uiPath = fullfile(launcherDir, 'phase_stamps_beta_ui', 'index.html');
            obj.HTML = uihtml(grid, 'HTMLSource', uiPath);
            obj.HTML.Layout.Row = 1;
            obj.HTML.Layout.Column = 1;
            obj.HTML.HTMLEventReceivedFcn = @(~,event) obj.onHtmlEvent(event);
            obj.configureLiveLog();
            obj.createTsPickerOverlay();
            obj.UIFigure.AutoResizeChildren = 'off';
            obj.UIFigure.SizeChangedFcn = @(~,~) obj.layoutTsPickerOverlay();

            for k = 1:numel(messages)
                obj.appendLog(messages{k});
            end
            if ~info.exists
                obj.appendLog(['No input_StaMPS.mat was found. PHASE loaded an initial ', ...
                    'configuration: review the detected values, select the StaMPS ', ...
                    'installation folder and press Save to create it.']);
            end
            if ~isempty(info.migratedFields)
                obj.appendLog(['Legacy input migrated in memory: ' ...
                    strjoin(info.migratedFields, ', ') '. Press Save to persist.']);
                obj.IsDirty = true;
            end
            obj.StatusDetail = ternary(obj.IsDirty, ...
                'Review detected values and save before starting', ...
                'Configuration loaded');
            drawnow;
            obj.sendState();
        end

        function appendLog(obj, message)
            message = char(string(message));
            timestamp = char(datetime('now', 'Format', 'HH:mm:ss'));
            entry = struct('time', timestamp, 'message', message);
            obj.Logs{end+1} = entry;
            fprintf('[PHASE beta %s] %s\n', timestamp, message);
            if ~isempty(obj.HTML) && isvalid(obj.HTML)
                try
                    sendEventToHTMLSource(obj.HTML, 'PhaseLog', entry);
                    drawnow limitrate
                catch
                end
            end
        end

        function openTsPicker(obj)
            try
                if isempty(obj.TSPickerOverlay) || ~isvalid(obj.TSPickerOverlay)
                    obj.createTsPickerOverlay();
                end
                obj.layoutTsPickerOverlay();
                obj.TSPickerOverlay.Visible = 'on';
                try, uistack(obj.TSPickerOverlay,'top'); catch, end
                drawnow;
                phase_stamps_beta.openTsPicker( ...
                    obj.WorkDir, obj.Config, obj.TSPickerContainer);
                obj.appendLog('TS Points picker loaded inside PHASE StaMPS.');
            catch ME
                try, obj.TSPickerOverlay.Visible = 'off'; catch, end
                obj.showError('TS picker unavailable', ME.message);
                obj.appendLog(['TS picker failed: ' ME.message]);
            end
        end

        function delete(obj)
            obj.endLiveDiary();
            try
                if ~isempty(obj.UIFigure) && isvalid(obj.UIFigure)
                    obj.UIFigure.CloseRequestFcn = [];
                    obj.UIFigure.SizeChangedFcn = [];
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
                name = char(string(event.HTMLEventName));
                payload = event.HTMLEventData;
            catch ME
                obj.appendLog(['Malformed HTML event: ' ME.message]);
                return
            end

            try
                switch lower(name)
                    case 'ready'
                        obj.sendState();
                    case 'load'
                        obj.loadFromDisk();
                    case 'save'
                        obj.saveFromPayload(payload);
                    case 'start'
                        obj.startFromPayload(payload);
                    case 'browse'
                        obj.browseFromPayload(payload);
                    case 'opentspicker'
                        obj.updateFromPayload(payload);
                        obj.openTsPicker();
                    case 'openworkdir'
                        openFolder(obj.WorkDir);
                    case 'openerrorlog'
                        errorLog = fullfile(obj.WorkDir, 'PHASE_StaMPS_error.log');
                        if exist(errorLog, 'file') == 2
                            openFile(errorLog);
                        else
                            obj.showError('Error log not found', ...
                                'No PHASE_StaMPS_error.log exists in this processing folder.');
                        end
                    otherwise
                        obj.appendLog(['Unknown interface event: ' name]);
                end
            catch ME
                obj.IsRunning = false;
                obj.Status = 'error';
                obj.StatusDetail = ME.message;
                obj.appendLog(['Interface action failed [' ME.identifier ']: ' ME.message]);
                obj.showError('PHASE StaMPS', ME.message);
                obj.sendState();
            end
        end

        function loadFromDisk(obj)
            if obj.IsRunning, return; end
            [loaded, info] = phase_stamps_beta.loadConfig(obj.WorkDir);
            [detected, fields, messages] = ...
                phase_stamps_beta.autoDetectConfig(loaded, obj.WorkDir);
            obj.SavedConfig = loaded;
            obj.Config = detected;
            obj.AutoDetectedFields = fields;
            obj.IsDirty = ~phase_stamps_beta.configsEqual(loaded, detected) || ...
                ~isempty(info.migratedFields) || ~info.exists;
            obj.Status = 'idle';
            obj.StatusDetail = ternary(obj.IsDirty, ...
                'Loaded; detected or migrated values require Save', ...
                'Configuration loaded from input_StaMPS.mat');
            if info.exists
                obj.appendLog(['Loaded configuration: ' info.path]);
            else
                obj.appendLog(['No input_StaMPS.mat was found. Initial values are ready; ', ...
                    'press Save to create ' info.path '.']);
            end
            for k = 1:numel(messages), obj.appendLog(messages{k}); end
            obj.sendState();
        end

        function saveFromPayload(obj, payload)
            if obj.IsRunning, return; end
            candidate = obj.configFromPayload(payload);
            [errors, warnings] = phase_stamps_beta.validateConfig(candidate, false);
            if ~isempty(errors)
                error('PHASE_StaMPS_beta:invalidConfiguration', '%s', strjoin(errors, newline));
            end
            phase_stamps_beta.saveConfig(obj.WorkDir, candidate);
            obj.Config = candidate;
            obj.SavedConfig = candidate;
            obj.IsDirty = false;
            obj.Status = 'saved';
            obj.StatusDetail = 'Configuration saved and ready to start';
            obj.appendLog(['Saved configuration: ' fullfile(obj.WorkDir, 'input_StaMPS.mat')]);
            for k = 1:numel(warnings), obj.appendLog(['Warning: ' warnings{k}]); end
            obj.sendState();
        end

        function startFromPayload(obj, payload)
            if obj.IsRunning, return; end
            candidate = obj.configFromPayload(payload);
            [errors, warnings] = phase_stamps_beta.validateConfig(candidate, true);
            if ~isempty(errors)
                error('PHASE_StaMPS_beta:invalidConfiguration', '%s', strjoin(errors, newline));
            end
            if isempty(obj.SavedConfig) || ...
                    ~phase_stamps_beta.configsEqual(candidate, obj.SavedConfig)
                error('PHASE_StaMPS_beta:unsavedConfiguration', ...
                    ['The visible configuration differs from input_StaMPS.mat. ' ...
                     'Press Save before Start.']);
            end

            runtimeMessages = phase_stamps_beta.prepareRuntime(candidate);
            obj.Config = candidate;
            obj.IsRunning = true;
            obj.Status = 'running';
            obj.StatusDetail = sprintf('Running StaMPS steps %s to %s', ...
                candidate.stamps_first_step, candidate.stamps_last_step);
            obj.Logs = {};
            obj.beginLiveDiary();
            diaryCleanup = onCleanup(@() obj.endLiveDiary()); %#ok<NASGU>
            obj.sendState();
            for k = 1:numel(warnings), obj.appendLog(['Warning: ' warnings{k}]); end
            for k = 1:numel(runtimeMessages), obj.appendLog(runtimeMessages{k}); end
            obj.appendLog(sprintf('Starting StaMPS steps %s -> %s in %s', ...
                candidate.stamps_first_step, candidate.stamps_last_step, obj.WorkDir));
            drawnow;

            adapter = phase_stamps_beta.LegacyAppAdapter(obj, candidate);
            try
                result = phase_stamps_beta.runProcessing(adapter);
            catch ME
                result = struct('ok', false, 'message', ME.message, ...
                    'identifier', ME.identifier);
                obj.appendLog(['Processing engine failed before its internal error handler [' ...
                    ME.identifier ']: ' ME.message]);
            end
            obj.IsRunning = false;
            if result.ok
                obj.Status = 'success';
                obj.StatusDetail = 'StaMPS processing completed';
            else
                obj.Status = 'error';
                obj.StatusDetail = result.message;
            end
            clear diaryCleanup
            obj.sendState();
            if result.ok
                obj.offerModelLaunch();
            end
        end

        function offerModelLaunch(obj)
            try
                choice = uiconfirm(obj.UIFigure, ...
                    ['StaMPS processing completed successfully. ' ...
                     'Do you want to open PHASE Model now?'], ...
                    'StaMPS completed', ...
                    'Options', {'Open PHASE Model','Not now'}, ...
                    'DefaultOption', 1, 'CancelOption', 2, 'Icon', 'success');
                if ~strcmp(choice, 'Open PHASE Model')
                    obj.appendLog('PHASE Model launch postponed.');
                    return
                end

                projectRoot = fileparts(obj.LauncherDir);
                modelLauncher = fullfile(projectRoot, 'PHASE_Model_beta.m');
                if ~isfile(modelLauncher)
                    error('PHASE_StaMPS_beta:modelLauncherMissing', ...
                        'PHASE Model launcher was not found: %s', modelLauncher);
                end
                addpath(projectRoot);
                obj.appendLog('Opening PHASE Model…');
                PHASE_Model_beta();
            catch ME
                obj.appendLog(['PHASE Model could not be opened [' ...
                    ME.identifier ']: ' ME.message]);
                obj.showError('PHASE Model unavailable', ME.message);
            end
        end

        function browseFromPayload(obj, payload)
            if obj.IsRunning, return; end
            obj.updateFromPayload(payload);
            if ~isstruct(payload) || ~isfield(payload, 'field')
                error('PHASE_StaMPS_beta:missingBrowseField', 'Browse action has no target field.');
            end
            fieldName = char(string(payload.field));
            if ~any(strcmp(fieldName, {'installation_folder','project_path'}))
                error('PHASE_StaMPS_beta:invalidBrowseField', 'Unsupported path field: %s', fieldName);
            end
            startFolder = obj.WorkDir;
            if isfield(obj.Config, fieldName) && isfolder(obj.Config.(fieldName))
                startFolder = obj.Config.(fieldName);
            end
            selected = uigetdir(startFolder, ['Select ' strrep(fieldName, '_', ' ')]);
            if ~isequal(selected, 0)
                obj.Config.(fieldName) = selected;
                obj.IsDirty = true;
                obj.Status = 'idle';
                obj.StatusDetail = 'Unsaved changes';
                obj.sendState();
            end
        end

        function updateFromPayload(obj, payload)
            obj.Config = obj.configFromPayload(payload);
            obj.IsDirty = isempty(obj.SavedConfig) || ...
                ~phase_stamps_beta.configsEqual(obj.Config, obj.SavedConfig);
        end

        function cfg = configFromPayload(obj, payload)
            if ~isstruct(payload) || ~isfield(payload, 'config')
                error('PHASE_StaMPS_beta:missingConfig', 'The interface did not send a configuration.');
            end
            cfg = phase_stamps_beta.uiToConfig(payload.config, obj.Config);
        end

        function sendState(obj)
            if isempty(obj.HTML) || ~isvalid(obj.HTML), return; end
            state = struct();
            state.kind = 'state';
            state.version = '6.0.0';
            state.workDir = obj.WorkDir;
            state.configPath = fullfile(obj.WorkDir, 'input_StaMPS.mat');
            state.schema = phase_stamps_beta.schema();
            state.config = phase_stamps_beta.configToUi(obj.Config);
            state.detectedFields = obj.AutoDetectedFields;
            state.dirty = obj.IsDirty;
            state.running = obj.IsRunning;
            state.status = obj.Status;
            state.statusDetail = obj.StatusDetail;
            state.logs = obj.Logs;
            state.liveLogUrl = obj.LiveLogUrl;
            obj.HTML.Data = state;
            drawnow limitrate
        end

        function configureLiveLog(obj)
            runtimeDir = fullfile(obj.LauncherDir, ...
                'phase_stamps_beta_ui','runtime_logs');
            if ~isfolder(runtimeDir), mkdir(runtimeDir); end
            token = char(java.util.UUID.randomUUID());
            fileName = ['stamps_' token '.log'];
            obj.LiveLogFile = fullfile(runtimeDir,fileName);
            obj.LiveLogUrl = ['runtime_logs/' fileName];
        end

        function beginLiveDiary(obj)
            obj.endLiveDiary();
            if isempty(obj.LiveLogFile), obj.configureLiveLog(); end
            try
                if isfile(obj.LiveLogFile), delete(obj.LiveLogFile); end
                diary(obj.LiveLogFile);
                diary on
                obj.DiaryActive = true;
            catch ME
                obj.DiaryActive = false;
                obj.appendLog(['Live Command Window capture unavailable: ' ME.message]);
            end
        end

        function endLiveDiary(obj)
            if ~obj.DiaryActive, return; end
            try, diary off; catch, end
            obj.DiaryActive = false;
        end

        function createTsPickerOverlay(obj)
            if ~isempty(obj.TSPickerOverlay) && isvalid(obj.TSPickerOverlay)
                return
            end
            obj.TSPickerOverlay = uipanel(obj.UIFigure, ...
                'BorderType','none','BackgroundColor',[0.985 0.988 0.994], ...
                'Visible','off');
            outer = uigridlayout(obj.TSPickerOverlay,[2 1]);
            outer.RowHeight = {58,'1x'};
            outer.Padding = [18 14 18 18];
            outer.RowSpacing = 10;

            header = uigridlayout(outer,[1 3]);
            header.Layout.Row = 1;
            header.ColumnWidth = {'1x','fit','fit'};
            header.Padding = [0 0 0 0];
            title = uilabel(header,'Text','TS Points', ...
                'FontName','Helvetica','FontSize',20,'FontWeight','bold', ...
                'FontColor',[0.27 0.275 0.275]);
            title.Layout.Column = 1;
            refresh = uibutton(header,'push','Text','Reload picker', ...
                'ButtonPushedFcn',@(~,~) obj.openTsPicker(), ...
                'BackgroundColor',[0.92 0.945 1.0], ...
                'FontColor',[0.208 0.396 0.812]);
            refresh.Layout.Column = 2;
            back = uibutton(header,'push','Text','Back to PHASE', ...
                'ButtonPushedFcn',@(~,~) obj.closeTsPicker(), ...
                'BackgroundColor',[1 1 1], ...
                'FontColor',[0.27 0.275 0.275]);
            back.Layout.Column = 3;

            obj.TSPickerContainer = uipanel(outer, ...
                'BorderType','line','BackgroundColor',[1 1 1]);
            obj.TSPickerContainer.Layout.Row = 2;
            obj.layoutTsPickerOverlay();
        end

        function layoutTsPickerOverlay(obj)
            if isempty(obj.TSPickerOverlay) || ~isvalid(obj.TSPickerOverlay) || ...
                    isempty(obj.UIFigure) || ~isvalid(obj.UIFigure)
                return
            end
            position = obj.UIFigure.Position;
            sidebarWidth = 254;
            if position(3) <= 1150, sidebarWidth = 224; end
            topbarHeight = 72;
            obj.TSPickerOverlay.Position = [sidebarWidth 0 ...
                max(100,position(3)-sidebarWidth) ...
                max(100,position(4)-topbarHeight)];
        end

        function closeTsPicker(obj)
            if ~isempty(obj.TSPickerOverlay) && isvalid(obj.TSPickerOverlay)
                obj.TSPickerOverlay.Visible = 'off';
            end
        end

        function showError(obj, titleText, message)
            try
                uialert(obj.UIFigure, char(string(message)), titleText);
            catch
                warning('%s: %s', titleText, char(string(message)));
            end
        end
    end
end

function pos = centeredPosition(width, height)
screen = get(groot, 'ScreenSize');
width = min(width, max(1000, screen(3) - 80));
height = min(height, max(700, screen(4) - 120));
pos = [max(20, (screen(3)-width)/2), max(40, (screen(4)-height)/2), width, height];
end

function out = ternary(condition, yesValue, noValue)
if condition, out = yesValue; else, out = noValue; end
end

function openFolder(pathValue)
if ispc
    winopen(pathValue);
else
    desktop = java.awt.Desktop.getDesktop();
    desktop.open(java.io.File(pathValue));
end
end

function openFile(pathValue)
if ispc
    winopen(pathValue);
else
    desktop = java.awt.Desktop.getDesktop();
    desktop.open(java.io.File(pathValue));
end
end
