classdef App < handle
    %APP One MATLAB window for PHASE projects and the three beta modules.

    properties (SetAccess = private)
        UIFigure
        TabGroup
        HomeTab
        PreprocessingTab
        StampsTab
        ModelTab
        ProjectLabel
        UpdateButton
        HomeTitle
        HomeDetails
        HomeStatus
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
            obj.ProjectLabel.Text = [char(string(project.name)) '  ·  ' paths.root];
            obj.updateHome();
            obj.refreshDatasets();
            obj.TabGroup.SelectedTab = obj.HomeTab;
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
            obj.TabGroup.SelectedTab = tab;
            obj.activateTab(tab);
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
            shell.RowHeight = {72,'1x'};
            shell.Padding = [0 0 0 0];
            shell.RowSpacing = 0;
            toolbar = uigridlayout(shell,[1 6]);
            toolbar.Layout.Row = 1;
            toolbar.ColumnWidth = {42,145,'1x',150,150,150};
            toolbar.Padding = [24 12 24 12];
            toolbar.ColumnSpacing = 10;
            toolbar.BackgroundColor = [1 1 1];
            logoPath = fullfile(obj.InstallRoot,'Logo_square.png');
            if isfile(logoPath)
                logo = uiimage(toolbar,'ImageSource',logoPath, ...
                    'ScaleMethod','fit');
                logo.Layout.Column = 1;
            end
            title = uilabel(toolbar,'Text','PHASE', ...
                'FontSize',25,'FontWeight','bold','FontColor',[53 101 207]/255);
            title.Layout.Column = 2;
            obj.ProjectLabel = uilabel(toolbar,'Text','No project open', ...
                'FontSize',13,'FontColor',[0.31 0.36 0.45]);
            obj.ProjectLabel.Layout.Column = 3;
            openButton = uibutton(toolbar,'push','Text','Open project', ...
                'ButtonPushedFcn',@(~,~) obj.chooseProject());
            openButton.Layout.Column = 4;
            styleButton(openButton,false);
            newButton = uibutton(toolbar,'push','Text','New project', ...
                'ButtonPushedFcn',@(~,~) obj.createProject());
            newButton.Layout.Column = 5;
            styleButton(newButton,true);
            obj.UpdateButton = uibutton(toolbar,'push','Text','Check for updates', ...
                'ButtonPushedFcn',@(~,~) obj.checkForUpdates());
            obj.UpdateButton.Layout.Column = 6;
            styleButton(obj.UpdateButton,false);
            if ~isfile(fullfile(fileparts(obj.InstallRoot),'install.json'))
                obj.UpdateButton.Enable = 'off';
                obj.UpdateButton.Tooltip = 'Available in managed PHASE installations.';
            end

            obj.TabGroup = uitabgroup(shell);
            obj.TabGroup.Layout.Row = 2;
            obj.HomeTab = uitab(obj.TabGroup,'Title','Project');
            obj.PreprocessingTab = uitab(obj.TabGroup,'Title','1 · Preprocessing');
            obj.StampsTab = uitab(obj.TabGroup,'Title','2 · StaMPS');
            obj.ModelTab = uitab(obj.TabGroup,'Title','3 · Model');
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
            obj.TabGroup.SelectionChangedFcn = @(~,event) obj.activateTab(event.NewValue);
        end

        function buildHome(obj)
            layout = uigridlayout(obj.HomeTab,[5 1]);
            layout.RowHeight = {132,26,238,118,'1x'};
            layout.Padding = [28 28 28 28];
            layout.RowSpacing = 16;
            layout.BackgroundColor = [0.975 0.978 0.985];

            welcome = uipanel(layout,'BorderType','none', ...
                'BackgroundColor',[0.92 0.945 1]);
            welcome.Layout.Row = 1;
            welcomeGrid = uigridlayout(welcome,[2 1]);
            welcomeGrid.RowHeight = {48,'1x'};
            welcomeGrid.Padding = [26 20 26 20];
            welcomeGrid.RowSpacing = 0;
            welcomeGrid.BackgroundColor = welcome.BackgroundColor;
            obj.HomeTitle = uilabel(welcomeGrid,'Text','Your PHASE workspace', ...
                'FontSize',28,'FontWeight','bold', ...
                'FontColor',[0.08 0.16 0.3]);
            obj.HomeTitle.Layout.Row = 1;
            obj.HomeDetails = uilabel(welcomeGrid, ...
                'Text','Create a project or open an existing one to get started.', ...
                'FontSize',15,'FontColor',[0.2 0.28 0.4]);
            obj.HomeDetails.Layout.Row = 2;

            workflow = uilabel(layout,'Text','WORKFLOW', ...
                'FontSize',11,'FontWeight','bold', ...
                'FontColor',[0.35 0.43 0.56]);
            workflow.Layout.Row = 2;
            cards = uigridlayout(layout,[1 3]);
            cards.Layout.Row = 3;
            cards.ColumnWidth = {'1x','1x','1x'};
            cards.ColumnSpacing = 14;
            cards.Padding = [0 0 0 0];
            cards.BackgroundColor = layout.BackgroundColor;
            obj.addWorkflowCard(cards,1,'01  PREPROCESSING', ...
                'Prepare SAR data', ...
                'Set your area of interest, process the image stack and export StaMPS inputs.', ...
                'Open Preprocessing',obj.PreprocessingTab);
            obj.addWorkflowCard(cards,2,'02  STAMPS', ...
                'Run PSI analysis', ...
                'Select a processed dataset, run StaMPS and export displacement results.', ...
                'Open StaMPS',obj.StampsTab);
            obj.addWorkflowCard(cards,3,'03  MODEL', ...
                'Explore results', ...
                'Build spatial and temporal models, figures and GIS-ready outputs.', ...
                'Open Model',obj.ModelTab);

            statusPanel = uipanel(layout,'BorderType','none', ...
                'BackgroundColor',[1 1 1]);
            statusPanel.Layout.Row = 4;
            statusGrid = uigridlayout(statusPanel,[2 1]);
            statusGrid.RowHeight = {28,'1x'};
            statusGrid.Padding = [22 17 22 17];
            statusGrid.RowSpacing = 2;
            statusGrid.BackgroundColor = [1 1 1];
            statusTitle = uilabel(statusGrid,'Text','PROJECT STATUS', ...
                'FontSize',11,'FontWeight','bold', ...
                'FontColor',[0.35 0.43 0.56]);
            statusTitle.Layout.Row = 1;
            obj.HomeStatus = uilabel(statusGrid, ...
                'Text','No project open. Your processing files and results will stay in the project folder.', ...
                'FontSize',14,'WordWrap','on', ...
                'FontColor',[0.2 0.28 0.4]);
            obj.HomeStatus.Layout.Row = 2;
        end

        function addWorkflowCard(obj,parent,column,step,heading,description,action,tab)
            card = uipanel(parent,'BorderType','line', ...
                'BackgroundColor',[1 1 1]);
            card.Layout.Column = column;
            grid = uigridlayout(card,[4 1]);
            grid.RowHeight = {23,37,'1x',42};
            grid.Padding = [22 20 22 18];
            grid.RowSpacing = 4;
            grid.BackgroundColor = [1 1 1];
            tag = uilabel(grid,'Text',step,'FontSize',11, ...
                'FontWeight','bold','FontColor',[53 101 207]/255);
            tag.Layout.Row = 1;
            headingLabel = uilabel(grid,'Text',heading,'FontSize',19, ...
                'FontWeight','bold','FontColor',[0.08 0.16 0.3]);
            headingLabel.Layout.Row = 2;
            descriptionLabel = uilabel(grid,'Text',description, ...
                'FontSize',13,'WordWrap','on', ...
                'FontColor',[0.31 0.36 0.45]);
            descriptionLabel.Layout.Row = 3;
            button = uibutton(grid,'push','Text',action, ...
                'ButtonPushedFcn',@(~,~) obj.goToTab(tab));
            button.Layout.Row = 4;
            styleButton(button,false);
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

        function activateTab(obj, tab)
            if isempty(obj.ProjectRoot) && tab ~= obj.HomeTab
                obj.TabGroup.SelectedTab = obj.HomeTab;
                uialert(obj.UIFigure,'Open or create a PHASE project first.', ...
                    'Project required');
                return
            end
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
            obj.TabGroup.SelectedTab = tab;
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
            if isempty(obj.HomeTitle) || ~isvalid(obj.HomeTitle), return; end
            if isempty(obj.ProjectRoot)
                obj.HomeTitle.Text = 'Your PHASE workspace';
                obj.HomeDetails.Text = 'Create a project or open an existing one to get started.';
                obj.HomeStatus.Text = [ ...
                    'No project open. Your processing files and results will stay in the project folder.'];
            else
                obj.HomeTitle.Text = char(string(obj.Project.name));
                p = phase_project.paths(obj.ProjectRoot);
                datasets = [dir(fullfile(p.stamps,'ASC_*')); ...
                    dir(fullfile(p.stamps,'DSC_*')); ...
                    dir(fullfile(p.stamps,'DES_*'))];
                count = nnz([datasets.isdir]);
                obj.HomeDetails.Text = 'One project for Preprocessing, StaMPS and Model.';
                obj.HomeStatus.Text = sprintf('Folder: %s    |    StaMPS datasets: %d', ...
                    obj.ProjectRoot,count);
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
