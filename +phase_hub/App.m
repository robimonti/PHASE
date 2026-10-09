classdef App < handle
    %APP One MATLAB window for PHASE projects and the three beta modules.

    properties (SetAccess = private)
        UIFigure
        Navigation
        HomeView
        HomeTab
        PreprocessingTab
        StampsTab
        ModelTab
        PreprocessingHost
        StampsHost
        ModelHost
        StampsSelector
        StampsPlaceholder
        PreprocessingApp = []
        StampsApp = []
        ModelApp = []
        InstallRoot
        ProjectRoot = ''
        Project = struct()
        CurrentStampsDir = ''
        ActiveSection = 'project'
        IsClosing = false
    end

    methods
        function obj = App(installRoot, projectRoot)
            if nargin < 1 || isempty(installRoot)
                installRoot = phase_project.installationRoot();
            end
            obj.InstallRoot = char(java.io.File(installRoot).getCanonicalPath());
            if nargin >= 2 && ~isempty(projectRoot)
                phase_project.open(projectRoot);
            end
            obj.buildWindow();
            if nargin >= 2 && ~isempty(projectRoot)
                obj.openProject(projectRoot);
            else
                obj.updateHome();
            end
        end

        function openProject(obj, projectRoot)
            [project, paths] = phase_project.open(projectRoot);
            if strcmp(paths.root,obj.ProjectRoot), return; end
            obj.assertCanSwitchProject();
            obj.unloadModules();
            obj.ProjectRoot = paths.root;
            obj.Project = project;
            obj.UIFigure.Name = ['PHASE · ' char(string(project.name))];
            obj.updateHome();
            obj.refreshDatasets();
            obj.showSection('project');
        end

        function refreshDatasets(obj)
            if isempty(obj.ProjectRoot)
                names = {'No project open'};
                paths = {''};
            else
                p = phase_project.paths(obj.ProjectRoot);
                datasets = [dir(fullfile(p.stamps,'ASC_*')); ...
                    dir(fullfile(p.stamps,'DSC_*')); ...
                    dir(fullfile(p.stamps,'DES_*'))];
                datasets = datasets([datasets.isdir]);
                names = {datasets.name};
                paths = cellfun(@(name) fullfile(p.stamps,name),names, ...
                    'UniformOutput',false);
                if isempty(names)
                    names = {'No dataset — complete Preprocessing first'};
                    paths = {''};
                end
            end
            previous = obj.CurrentStampsDir;
            obj.StampsSelector.Items = names;
            obj.StampsSelector.ItemsData = paths;
            if ~isempty(previous) && any(strcmp(paths,previous))
                obj.StampsSelector.Value = previous;
            else
                obj.StampsSelector.Value = paths{1};
            end
            obj.StampsPlaceholder.Text = obj.stampsHelpText();
            obj.StampsSelector.Enable = ternary(~isempty(paths{1}),'on','off');
            obj.updateHome();
        end

        function showSection(obj, name)
            switch lower(char(string(name)))
                case {'home','project','progetto'}
                    tab = obj.HomeTab;
                case {'preprocessing','preprocess'}
                    tab = obj.PreprocessingTab;
                case {'stamps','psi'}
                    tab = obj.StampsTab;
                case {'model','modello'}
                    tab = obj.ModelTab;
                otherwise
                    error('PHASE:HubSectionUnknown', ...
                        'Unknown PHASE section: %s',char(string(name)));
            end
            obj.selectSection(tab);
        end

        function delete(obj)
            if obj.IsClosing, return; end
            obj.IsClosing = true;
            try
                if ~isempty(obj.UIFigure) && isvalid(obj.UIFigure)
                    obj.UIFigure.CloseRequestFcn = [];
                end
                obj.unloadModules();
                if ~isempty(obj.UIFigure) && isvalid(obj.UIFigure)
                    obj.UIFigure.UserData = [];
                    delete(obj.UIFigure);
                end
            catch
            end
        end
    end

    methods (Access = private)
        function buildWindow(obj)
            obj.UIFigure = uifigure('Name','PHASE · Workspace', ...
                'Color',[1 1 1], ...
                'Position',centeredPosition(1540,960));
            obj.UIFigure.UserData = obj;
            obj.UIFigure.CloseRequestFcn = @(~,~) obj.requestClose();

            shell = uigridlayout(obj.UIFigure,[2 1]);
            shell.RowHeight = {130,'1x'};
            shell.Padding = [0 0 0 0];
            shell.RowSpacing = 0;
            htmlSource = fullfile(obj.InstallRoot,'PHASE_Hub_UI.html');
            if ~isfile(htmlSource)
                error('PHASE:HubViewMissing','The PHASE hub interface is missing: %s',htmlSource);
            end
            obj.Navigation = uihtml(shell,'HTMLSource',htmlSource, ...
                'DataChangedFcn',@(src,event) obj.handleHtmlAction(src,event));
            obj.Navigation.Layout.Row = 1;
            obj.Navigation.Data = struct('view','nav','active','project', ...
                'version',obj.displayVersion(), ...
                'updateEnabled',isfile(fullfile(fileparts(obj.InstallRoot),'install.json')));

            obj.HomeTab = uipanel(shell,'BorderType','none','BackgroundColor',[1 1 1]);
            obj.PreprocessingTab = uipanel(shell,'BorderType','none','BackgroundColor',[1 1 1]);
            obj.StampsTab = uipanel(shell,'BorderType','none','BackgroundColor',[1 1 1]);
            obj.ModelTab = uipanel(shell,'BorderType','none','BackgroundColor',[1 1 1]);
            panels = {obj.HomeTab,obj.PreprocessingTab,obj.StampsTab,obj.ModelTab};
            for k = 1:numel(panels)
                panels{k}.Layout.Row = 2;
                panels{k}.Layout.Column = 1;
                panels{k}.Visible = ternary(k == 1,'on','off');
            end
            obj.buildHome();
            obj.PreprocessingHost = uipanel(obj.PreprocessingTab, ...
                'BorderType','none','BackgroundColor',[1 1 1], ...
                'Position',[0 0 1480 830]);
            obj.fillTab(obj.PreprocessingTab,obj.PreprocessingHost);
            obj.buildStampsTab();
            obj.ModelHost = uipanel(obj.ModelTab, ...
                'BorderType','none','BackgroundColor',[1 1 1], ...
                'Position',[0 0 1480 830]);
            obj.fillTab(obj.ModelTab,obj.ModelHost);
        end

        function buildHome(obj)
            layout = uigridlayout(obj.HomeTab,[1 1]);
            layout.Padding = [0 0 0 0];
            obj.HomeView = uihtml(layout, ...
                'HTMLSource',fullfile(obj.InstallRoot,'PHASE_Hub_UI.html'), ...
                'DataChangedFcn',@(src,event) obj.handleHtmlAction(src,event));
            obj.HomeView.Layout.Row = 1;
            obj.HomeView.Layout.Column = 1;
        end

        function buildStampsTab(obj)
            layout = uigridlayout(obj.StampsTab,[2 1]);
            layout.RowHeight = {48,'1x'};
            layout.Padding = [0 0 0 0];
            layout.RowSpacing = 0;
            bar = uigridlayout(layout,[1 3]);
            bar.Layout.Row = 1;
            bar.ColumnWidth = {140,'1x',115};
            bar.Padding = [15 6 15 6];
            bar.BackgroundColor = [0.94 0.96 0.99];
            label = uilabel(bar,'Text','Dataset StaMPS', ...
                'FontWeight','bold');
            label.Layout.Column = 1;
            obj.StampsSelector = uidropdown(bar, ...
                'Items',{'No project open'},'ItemsData',{''}, ...
                'ValueChangedFcn',@(~,~) obj.loadStamps());
            obj.StampsSelector.Layout.Column = 2;
            refresh = uibutton(bar,'push','Text','Refresh', ...
                'ButtonPushedFcn',@(~,~) obj.refreshDatasets());
            refresh.Layout.Column = 3;
            styleButton(refresh,false);
            obj.StampsHost = uipanel(layout,'BorderType','none', ...
                'BackgroundColor',[1 1 1]);
            obj.StampsHost.Layout.Row = 2;
            obj.StampsPlaceholder = uilabel(obj.StampsHost, ...
                'Text','Open a project and select a StaMPS dataset.', ...
                'FontSize',17,'HorizontalAlignment','center', ...
                'FontColor',[0.28 0.34 0.45], ...
                'Position',[260 330 900 90]);
        end

        function fillTab(~,tab,panel)
            layout = uigridlayout(tab,[1 1]);
            layout.Padding = [0 0 0 0];
            panel.Parent = layout;
            panel.Layout.Row = 1;
            panel.Layout.Column = 1;
        end

        function handleHtmlAction(obj,~,event)
            data = event.Data;
            if ~isstruct(data) || ~isfield(data,'action'), return; end
            action = char(string(data.action));
            switch action
                case 'open'
                    obj.chooseProject();
                case 'new'
                    obj.createProject();
                case 'updates'
                    obj.checkForUpdates();
                case {'project','preprocessing','stamps','model'}
                    obj.showSection(action);
            end
        end

        function selectSection(obj,tab)
            if isempty(obj.ProjectRoot) && tab ~= obj.HomeTab
                uialert(obj.UIFigure,'Open or create a PHASE project first.', ...
                    'Project required');
                return
            end
            panels = {obj.HomeTab,obj.PreprocessingTab,obj.StampsTab,obj.ModelTab};
            names = {'project','preprocessing','stamps','model'};
            for k = 1:numel(panels)
                panels{k}.Visible = ternary(panels{k} == tab,'on','off');
                if panels{k} == tab, obj.ActiveSection = names{k}; end
            end
            obj.Navigation.Data = struct('view','nav', ...
                'active',obj.ActiveSection,'version',obj.displayVersion(), ...
                'updateEnabled',isfile(fullfile(fileparts(obj.InstallRoot),'install.json')));
            obj.activateTab(tab);
        end

        function value = displayVersion(obj)
            value = 'v7.0.0 preview';
            metadata = fullfile(fileparts(obj.InstallRoot),'install.json');
            if ~isfile(metadata), return; end
            try
                info = jsondecode(fileread(metadata));
                if isfield(info,'version')
                    installed = char(string(info.version));
                    if ~isempty(installed) && ~strcmp(installed,'dev')
                        if startsWith(installed,'v')
                            value = installed;
                        else
                            value = ['v' installed];
                        end
                    end
                end
            catch
            end
        end

        function activateTab(obj, tab)
            try
                if tab == obj.PreprocessingTab
                    if isempty(obj.PreprocessingApp) || ~isvalid(obj.PreprocessingApp)
                        obj.PreprocessingApp = phase_preprocessing_beta.App( ...
                            obj.InstallRoot,obj.ProjectRoot,obj.PreprocessingHost);
                    end
                elseif tab == obj.ModelTab
                    if isempty(obj.ModelApp) || ~isvalid(obj.ModelApp)
                        obj.ModelApp = phase_model_beta.App( ...
                            obj.InstallRoot,obj.ProjectRoot,obj.ModelHost);
                    end
                elseif tab == obj.StampsTab
                    obj.refreshDatasets();
                    obj.loadStamps();
                else
                    obj.updateHome();
                end
            catch ME
                uialert(obj.UIFigure,ME.message,'Cannot open section');
            end
        end

        function loadStamps(obj)
            if isempty(obj.ProjectRoot), return; end
            selected = char(string(obj.StampsSelector.Value));
            if isempty(selected)
                obj.StampsPlaceholder.Visible = 'on';
                return
            end
            if strcmp(selected,obj.CurrentStampsDir) && ...
                    ~isempty(obj.StampsApp) && isvalid(obj.StampsApp)
                return
            end
            if ~isempty(obj.StampsApp) && isvalid(obj.StampsApp)
                if obj.StampsApp.IsRunning
                    obj.StampsSelector.Value = obj.CurrentStampsDir;
                    uialert(obj.UIFigure, ...
                        'Wait for StaMPS processing to finish before switching datasets.', ...
                        'Processing in progress');
                    return
                end
                delete(obj.StampsApp);
                obj.StampsApp = [];
            end
            obj.StampsPlaceholder.Visible = 'off';
            try
                obj.StampsApp = phase_stamps_beta.App(selected, ...
                    fullfile(obj.InstallRoot,'PHASE_Preprocessing'),obj.StampsHost);
                obj.CurrentStampsDir = selected;
            catch ME
                obj.CurrentStampsDir = '';
                obj.StampsPlaceholder.Text = ['Cannot open dataset: ' ME.message];
                obj.StampsPlaceholder.Visible = 'on';
                rethrow(ME)
            end
        end

        function goToTab(obj,tab)
            obj.selectSection(tab);
        end

        function chooseProject(obj)
            selected = uigetdir(pwd,'Select a PHASE project');
            if isequal(selected,0), return; end
            try
                obj.openProject(selected);
            catch ME
                uialert(obj.UIFigure,ME.message,'Cannot open project');
            end
        end

        function createProject(obj)
            selected = uigetdir(pwd,'Select an empty folder for the new project');
            if isequal(selected,0), return; end
            try
                obj.assertCanSwitchProject();
                [~,paths] = phase_project.create(selected);
                obj.openProject(paths.root);
            catch ME
                uialert(obj.UIFigure,ME.message,'Cannot create project');
            end
        end

        function checkForUpdates(obj)
            try
                result = obj.runUpdater('check');
                if ~logical(result.updateAvailable)
                    uialert(obj.UIFigure, ...
                        sprintf('Installed version: %s. No updates available.', ...
                        char(string(result.current))), 'PHASE is up to date');
                    return
                end
                if obj.hasActiveWork()
                    uialert(obj.UIFigure, ...
                        'Finish current processing before preparing an update.', ...
                        'Processing in progress');
                    return
                end
                choice = uiconfirm(obj.UIFigure, ...
                    sprintf('PHASE %s is available. Download it now? It will be installed at the next launch.', ...
                    char(string(result.available))), ...
                    'PHASE update', ...
                    'Options',{'Download','Cancel'}, ...
                    'DefaultOption','Download','CancelOption','Cancel');
                if ~strcmp(choice,'Download'), return; end
                progress = uiprogressdlg(obj.UIFigure, ...
                    'Title','PHASE update', ...
                    'Message','Downloading and verifying the release...', ...
                    'Indeterminate','on');
                try
                    prepared = obj.runUpdater('prepare');
                    delete(progress);
                catch ME
                    delete(progress);
                    rethrow(ME)
                end
                if logical(prepared.prepared)
                    uialert(obj.UIFigure, ...
                        ['Update ready. Close MATLAB and relaunch PHASE ' ...
                        'from the application shortcut to install it.'], ...
                        'Restart required');
                end
            catch ME
                uialert(obj.UIFigure,ME.message,'Update unavailable');
            end
        end

        function result = runUpdater(obj,action)
            prefix = fileparts(obj.InstallRoot);
            python = getenv('PHASE_PYTHON');
            if isempty(python) && ispc
                config = fullfile(getenv('APPDATA'),'PHASE','python.txt');
                if isfile(config), python = strtrim(fileread(config)); end
            end
            if isempty(python)
                if ispc, python = 'python'; else, python = 'python3'; end
            end
            script = fullfile(obj.InstallRoot,'phase_update.py');
            if ~isfile(script)
                error('PHASE:UpdaterMissing','PHASE updater was not found in this installation.');
            end
            command = sprintf('"%s" "%s" %s --prefix "%s"', ...
                python,script,action,prefix);
            [status,output] = system(command);
            try
                result = jsondecode(strtrim(output));
            catch
                error('PHASE:UpdaterResponse','Unexpected updater response: %s',output);
            end
            if status ~= 0
                if isfield(result,'error')
                    error('PHASE:UpdaterFailed','%s',char(string(result.error)));
                end
                error('PHASE:UpdaterFailed','Update failed.');
            end
        end

        function updateHome(obj)
            if isempty(obj.HomeView) || ~isvalid(obj.HomeView), return; end
            state = struct('view','home','projectName','', ...
                'projectPath','','inputCount','—', ...
                'datasetCount','—','exportCount','—');
            if isempty(obj.ProjectRoot)
                obj.HomeView.Data = state;
            else
                p = phase_project.paths(obj.ProjectRoot);
                datasets = [dir(fullfile(p.stamps,'ASC_*')); ...
                    dir(fullfile(p.stamps,'DSC_*')); ...
                    dir(fullfile(p.stamps,'DES_*'))];
                state.projectName = char(string(obj.Project.name));
                state.projectPath = obj.ProjectRoot;
                state.inputCount = countFiles(p.raw);
                state.datasetCount = nnz([datasets.isdir]);
                state.exportCount = countFiles(p.exports);
                obj.HomeView.Data = state;
            end
        end

        function textValue = stampsHelpText(obj)
            if isempty(obj.ProjectRoot)
                textValue = 'Open or create a PHASE project.';
            else
                textValue = ['No StaMPS dataset is available yet. Complete ' ...
                    'Preprocessing, then press Refresh.'];
            end
        end

        function assertCanSwitchProject(obj)
            if obj.hasActiveWork()
                error('PHASE:HubBusy', ...
                    'Processing is in progress. Finish or stop it before switching projects.');
            end
        end

        function active = hasActiveWork(obj)
            active = ~isempty(obj.PreprocessingApp) && ...
                isvalid(obj.PreprocessingApp) && obj.PreprocessingApp.IsRunning;
            active = active || (~isempty(obj.PreprocessingApp) && ...
                isvalid(obj.PreprocessingApp) && ...
                isstruct(obj.PreprocessingApp.Transfer) && ...
                obj.PreprocessingApp.Transfer.active);
            active = active || (~isempty(obj.StampsApp) && ...
                isvalid(obj.StampsApp) && obj.StampsApp.IsRunning);
            active = active || (~isempty(obj.ModelApp) && ...
                isvalid(obj.ModelApp) && obj.ModelApp.IsRunning);
        end

        function unloadModules(obj)
            if ~isempty(obj.StampsApp) && isvalid(obj.StampsApp)
                delete(obj.StampsApp);
            end
            obj.StampsApp = [];
            obj.CurrentStampsDir = '';
            if ~isempty(obj.PreprocessingApp) && isvalid(obj.PreprocessingApp)
                delete(obj.PreprocessingApp);
            end
            obj.PreprocessingApp = [];
            if ~isempty(obj.ModelApp) && isvalid(obj.ModelApp)
                delete(obj.ModelApp);
            end
            obj.ModelApp = [];
            if ~isempty(obj.StampsPlaceholder) && isvalid(obj.StampsPlaceholder)
                obj.StampsPlaceholder.Visible = 'on';
            end
        end

        function requestClose(obj)
            if obj.hasActiveWork()
                choice = uiconfirm(obj.UIFigure, ...
                    'Processing is still running. Close PHASE and stop it?', ...
                    'Processing in progress','Options',{'Cancel','Close PHASE'}, ...
                    'DefaultOption','Cancel','CancelOption','Cancel');
                if ~strcmp(choice,'Close PHASE'), return; end
            end
            delete(obj);
        end
    end
end

function position = centeredPosition(width,height)
screen = get(groot,'ScreenSize');
width = min(width,max(1100,screen(3)-50));
height = min(height,max(760,screen(4)-80));
position = [max(20,(screen(3)-width)/2), ...
    max(30,(screen(4)-height)/2),width,height];
end

function value = ternary(condition,yesValue,noValue)
if condition, value = yesValue; else, value = noValue; end
end

function count = countFiles(folder)
count = 0;
if ~isfolder(folder), return; end
entries = dir(folder);
count = nnz(~[entries.isdir]);
end

function styleButton(button,isPrimary)
button.FontSize = 13;
button.FontWeight = 'bold';
if ispc, button.FontName = 'Segoe UI'; end
if isPrimary
    button.BackgroundColor = [53 101 207]/255;
    button.FontColor = [1 1 1];
else
    button.BackgroundColor = [0.92 0.945 1];
    button.FontColor = [53 101 207]/255;
end
end
