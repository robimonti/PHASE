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
                names = {'Nessun progetto'};
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
                    names = {'Nessun dataset — completa il Preprocessing'};
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
                'Color',[0.965 0.973 0.99], ...
                'Position',centeredPosition(1540,960));
            obj.UIFigure.UserData = obj;
            obj.UIFigure.CloseRequestFcn = @(~,~) obj.requestClose();

            shell = uigridlayout(obj.UIFigure,[2 1]);
            shell.RowHeight = {60,'1x'};
            shell.Padding = [0 0 0 0];
            shell.RowSpacing = 0;
            toolbar = uigridlayout(shell,[1 5]);
            toolbar.Layout.Row = 1;
            toolbar.ColumnWidth = {180,'1x',138,138,138};
            toolbar.Padding = [18 10 18 10];
            toolbar.BackgroundColor = [0.09 0.15 0.29];
            title = uilabel(toolbar,'Text','PHASE', ...
                'FontSize',25,'FontWeight','bold','FontColor',[1 1 1]);
            title.Layout.Column = 1;
            obj.ProjectLabel = uilabel(toolbar,'Text','Nessun progetto aperto', ...
                'FontSize',13,'FontColor',[0.84 0.9 1]);
            obj.ProjectLabel.Layout.Column = 2;
            openButton = uibutton(toolbar,'push','Text','Apri progetto', ...
                'ButtonPushedFcn',@(~,~) obj.chooseProject());
            openButton.Layout.Column = 3;
            newButton = uibutton(toolbar,'push','Text','Nuovo progetto', ...
                'ButtonPushedFcn',@(~,~) obj.createProject());
            newButton.Layout.Column = 4;
            obj.UpdateButton = uibutton(toolbar,'push','Text','Cerca update', ...
                'ButtonPushedFcn',@(~,~) obj.checkForUpdates());
            obj.UpdateButton.Layout.Column = 5;
            if ~isfile(fullfile(fileparts(obj.InstallRoot),'install.json'))
                obj.UpdateButton.Enable = 'off';
                obj.UpdateButton.Tooltip = 'Disponibile nelle installazioni PHASE 7.';
            end

            obj.TabGroup = uitabgroup(shell);
            obj.TabGroup.Layout.Row = 2;
            obj.HomeTab = uitab(obj.TabGroup,'Title','Progetto');
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
            layout = uigridlayout(obj.HomeTab,[4 1]);
            layout.RowHeight = {65,105,80,'1x'};
            layout.Padding = [38 36 38 36];
            layout.RowSpacing = 15;
            layout.BackgroundColor = [0.965 0.973 0.99];
            obj.HomeTitle = uilabel(layout,'Text','Il tuo workspace PHASE', ...
                'FontSize',29,'FontWeight','bold','FontColor',[0.08 0.16 0.3]);
            obj.HomeTitle.Layout.Row = 1;
            obj.HomeDetails = uilabel(layout, ...
                'Text','Apri o crea un progetto per iniziare.', ...
                'FontSize',15,'FontColor',[0.2 0.28 0.4]);
            obj.HomeDetails.Layout.Row = 2;
            actions = uigridlayout(layout,[1 3]);
            actions.Layout.Row = 3;
            actions.ColumnWidth = {'1x','1x','1x'};
            actions.Padding = [0 0 0 0];
            names = {'1 · Preprocessing','2 · StaMPS','3 · Model'};
            tabs = {obj.PreprocessingTab,obj.StampsTab,obj.ModelTab};
            for k = 1:3
                target = tabs{k};
                button = uibutton(actions,'push','Text',names{k}, ...
                    'FontSize',14,'ButtonPushedFcn',@(~,~) obj.goToTab(target));
                button.Layout.Column = k;
            end
            note = uilabel(layout, ...
                'Text',['Un solo progetto per tutte le sezioni. I file di lavoro e i ' ...
                'risultati restano nel progetto; il codice resta nell’installazione.'], ...
                'FontSize',13,'FontColor',[0.37 0.43 0.52]);
            note.Layout.Row = 4;
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
                'Items',{'Nessun progetto'},'ItemsData',{''}, ...
                'ValueChangedFcn',@(~,~) obj.loadStamps());
            obj.StampsSelector.Layout.Column = 2;
            refresh = uibutton(bar,'push','Text','Aggiorna', ...
                'ButtonPushedFcn',@(~,~) obj.refreshDatasets());
            refresh.Layout.Column = 3;
            obj.StampsHost = uipanel(layout,'BorderType','none', ...
                'BackgroundColor',[1 1 1]);
            obj.StampsHost.Layout.Row = 2;
            obj.StampsPlaceholder = uilabel(obj.StampsHost, ...
                'Text','Apri un progetto e scegli un dataset StaMPS.', ...
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
                uialert(obj.UIFigure,'Apri o crea prima un progetto PHASE.', ...
                    'Progetto richiesto');
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
                uialert(obj.UIFigure,ME.message,'Impossibile aprire la sezione');
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
                        'Attendi che l’elaborazione StaMPS termini prima di cambiare dataset.', ...
                        'Elaborazione in corso');
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
                obj.StampsPlaceholder.Text = ['Dataset non apribile: ' ME.message];
                obj.StampsPlaceholder.Visible = 'on';
                rethrow(ME)
            end
        end

        function goToTab(obj,tab)
            obj.TabGroup.SelectedTab = tab;
            obj.activateTab(tab);
        end

        function chooseProject(obj)
            selected = uigetdir(pwd,'Seleziona un progetto PHASE');
            if isequal(selected,0), return; end
            try
                obj.openProject(selected);
            catch ME
                uialert(obj.UIFigure,ME.message,'Impossibile aprire il progetto');
            end
        end

        function createProject(obj)
            selected = uigetdir(pwd,'Seleziona una cartella vuota per il progetto');
            if isequal(selected,0), return; end
            try
                obj.assertCanSwitchProject();
                [~,paths] = phase_project.create(selected);
                obj.openProject(paths.root);
            catch ME
                uialert(obj.UIFigure,ME.message,'Impossibile creare il progetto');
            end
        end

        function checkForUpdates(obj)
            try
                result = obj.runUpdater('check');
                if ~logical(result.updateAvailable)
                    uialert(obj.UIFigure, ...
                        sprintf('Versione installata: %s. Nessun aggiornamento disponibile.', ...
                        char(string(result.current))), 'PHASE aggiornato');
                    return
                end
                if obj.hasActiveWork()
                    uialert(obj.UIFigure, ...
                        'Termina le elaborazioni prima di preparare un aggiornamento.', ...
                        'Elaborazione in corso');
                    return
                end
                choice = uiconfirm(obj.UIFigure, ...
                    sprintf('È disponibile PHASE %s. Scaricarlo ora? Verrà installato al prossimo avvio.', ...
                    char(string(result.available))), ...
                    'Aggiornamento PHASE', ...
                    'Options',{'Scarica','Annulla'}, ...
                    'DefaultOption','Scarica','CancelOption','Annulla');
                if ~strcmp(choice,'Scarica'), return; end
                progress = uiprogressdlg(obj.UIFigure, ...
                    'Title','Aggiornamento PHASE', ...
                    'Message','Download e verifica della release...', ...
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
                        ['Aggiornamento pronto. Chiudi MATLAB e riapri PHASE ' ...
                        'dal collegamento dell’applicazione per applicarlo.'], ...
                        'Riavvio richiesto');
                end
            catch ME
                uialert(obj.UIFigure,ME.message,'Aggiornamento non disponibile');
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
                error('PHASE:UpdaterMissing','Updater PHASE non trovato nell’installazione.');
            end
            command = sprintf('"%s" "%s" %s --prefix "%s"', ...
                python,script,action,prefix);
            [status,output] = system(command);
            try
                result = jsondecode(strtrim(output));
            catch
                error('PHASE:UpdaterResponse','Risposta inattesa dell’updater: %s',output);
            end
            if status ~= 0
                if isfield(result,'error')
                    error('PHASE:UpdaterFailed','%s',char(string(result.error)));
                end
                error('PHASE:UpdaterFailed','Aggiornamento non riuscito.');
            end
        end

        function updateHome(obj)
            if isempty(obj.HomeTitle) || ~isvalid(obj.HomeTitle), return; end
            if isempty(obj.ProjectRoot)
                obj.HomeTitle.Text = 'Il tuo workspace PHASE';
                obj.HomeDetails.Text = 'Apri o crea un progetto per iniziare.';
            else
                obj.HomeTitle.Text = char(string(obj.Project.name));
                p = phase_project.paths(obj.ProjectRoot);
                datasets = [dir(fullfile(p.stamps,'ASC_*')); ...
                    dir(fullfile(p.stamps,'DSC_*')); ...
                    dir(fullfile(p.stamps,'DES_*'))];
                count = nnz([datasets.isdir]);
                obj.HomeDetails.Text = sprintf('Cartella: %s\nDataset StaMPS: %d', ...
                    obj.ProjectRoot,count);
            end
        end

        function textValue = stampsHelpText(obj)
            if isempty(obj.ProjectRoot)
                textValue = 'Apri o crea un progetto PHASE.';
            else
                textValue = ['Nessun dataset StaMPS disponibile. Completa il ' ...
                    'Preprocessing e premi Aggiorna.'];
            end
        end

        function assertCanSwitchProject(obj)
            if obj.hasActiveWork()
                error('PHASE:HubBusy', ...
                    'Un’elaborazione è in corso. Termina o arresta il lavoro prima di cambiare progetto.');
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
                    'Un’elaborazione è ancora in corso. Chiudere PHASE e interromperla?', ...
                    'Lavoro in corso','Options',{'Annulla','Chiudi PHASE'}, ...
                    'DefaultOption','Annulla','CancelOption','Annulla');
                if ~strcmp(choice,'Chiudi PHASE'), return; end
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
