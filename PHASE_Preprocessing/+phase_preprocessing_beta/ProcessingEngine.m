classdef ProcessingEngine < handle
    %PROCESSINGENGINE .m-only preprocessing runner without App Designer UI.
    % The visible controller saves configuration and owns the download map.

    properties
        ProjectRoot = ''
        InstallRoot = ''
        StopFlag = false
        ActiveProcess = []
        GenerateCoherence = true
        GenerateLia = true
        constellation = "SEN"
        update_processed_data_SEN = 0
        ExternalLogCallback = []
        ExternalProgressCallback = []
    end

    methods
        function app = ProcessingEngine(projectRoot,installRoot)
            if nargin >= 1, app.ProjectRoot = char(string(projectRoot)); end
            if nargin >= 2, app.InstallRoot = char(string(installRoot)); end
        end

        function updateOutput(app,message)
            message = char(string(message));
            if ~isempty(app.ExternalLogCallback)
                app.ExternalLogCallback(message);
            else
                fprintf('%s\n',message);
            end
        end

        function notifyBetaProgress(app,percentage,phase)
            if isempty(app.ExternalProgressCallback), return; end
            app.ExternalProgressCallback(struct( ...
                'percentage',double(percentage), ...
                'phase',char(string(phase)),'indeterminate',false));
        end

        function stopped = forceStopActiveProcess(app)
            stopped = false;
            if isempty(app.ActiveProcess), return; end
            try
                if app.ActiveProcess.isAlive()
                    stopped = phase_preprocessing_beta.forceStopProcess(app.ActiveProcess);
                end
            catch
            end
        end

        function deleteFileIfExists(~,filePath)
            if isfile(filePath), delete(filePath); end
        end

        function reprojectGeoTiffIfGeographic(app,filePath,epsgCode,bandIndex)
            if nargin < 4, bandIndex = []; end
            normalized = lower(strrep(filePath,'\','/'));
            if (~app.GenerateCoherence && contains(normalized,'/coherence/')) || ...
                    (~app.GenerateLia && contains(normalized,'/lia/'))
                return
            end
            app.notifyBetaProgress(94,['Reprojecting ' shortFileName(filePath)]);
            phase_preprocessing_beta.reprojectGeoTiff(filePath,epsgCode,bandIndex);
            app.notifyBetaProgress(96,['Reprojected ' shortFileName(filePath)]);
        end

        function dateText = parseSENDateFromFilename(app, filename)
            dateText = '';
            token = regexp(filename, '_(\d{8})T\d{6}_', 'tokens', 'once');
            if ~isempty(token)
                dateText = token{1};
            end
        end

        function dateList = getSENRootZipDates(app, projectFolder)
            dateList = {};
            slavesFolder = fullfile(projectFolder, 'slaves');
            rootZips = dir(fullfile(slavesFolder, '*.zip'));

            for k = 1:numel(rootZips)
                dateText = app.parseSENDateFromFilename(rootZips(k).name);
                if ~isempty(dateText)
                    dateList{end + 1} = dateText; %#ok<AGROW>
                end
            end

            dateDirs = dir(slavesFolder);
            for k = 1:numel(dateDirs)
                dateText = dateDirs(k).name;
                if ~dateDirs(k).isdir || isempty(regexp(dateText, '^\d{8}$', 'once'))
                    continue
                end

                folderZips = dir(fullfile(dateDirs(k).folder, dateText, '*.zip'));
                if ~isempty(folderZips) && ~app.hasSENProcessedSlaveDate(projectFolder, dateText)
                    dateList{end + 1} = dateText; %#ok<AGROW>
                end
            end

            dateList = unique(dateList);
        end

        function isProcessed = hasSENProcessedSlaveDate(app, projectFolder, dateText)
            % Coregistration and interferogram products are the non-optional
            % slave outputs used to mark a date as already processed.
            coregFiles = dir(fullfile(projectFolder, 'coreg', ['*_' dateText '.dim']));
            ifgFiles = dir(fullfile(projectFolder, 'ifg', ['*_' dateText '.dim']));
            isProcessed = ~isempty(coregFiles) && ~isempty(ifgFiles);
        end

        function hiddenPaths = isolateSENUpdateWorkspace(app, projectFolder, updateDates)
            hiddenPaths = struct('original', {}, 'hidden', {});
            timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
            hiddenRoot = fullfile(projectFolder, '_update_hidden', timestamp);

            try
                slavesFolder = fullfile(projectFolder, 'slaves');
                hiddenPaths = app.hideDateDirectoriesExcept(hiddenPaths, slavesFolder, hiddenRoot, 'slaves', updateDates);

                splitFolder = fullfile(projectFolder, 'split');
                hiddenPaths = app.hideDateDirectoriesExcept(hiddenPaths, splitFolder, hiddenRoot, 'split', updateDates);

                productFolders = {'coreg', 'ifg'};
                for f = 1:numel(productFolders)
                    folderName = productFolders{f};
                    folderPath = fullfile(projectFolder, folderName);
                    if exist(folderPath, 'dir') ~= 7
                        continue
                    end

                    items = dir(folderPath);
                    for k = 1:numel(items)
                        name = items(k).name;
                        if strcmp(name, '.') || strcmp(name, '..')
                            continue
                        end

                        if ~app.nameContainsAnySENDate(name, updateDates)
                            originalPath = fullfile(items(k).folder, name);
                            relativePath = fullfile(folderName, name);
                            hiddenPaths = app.hideSENUpdatePath(hiddenPaths, originalPath, hiddenRoot, relativePath);
                        end
                    end
                end
            catch ME
                app.restoreSENUpdateWorkspace(hiddenPaths);
                rethrow(ME);
            end
        end

        function hiddenPaths = hideDateDirectoriesExcept(app, hiddenPaths, folderPath, hiddenRoot, groupName, updateDates)
            if exist(folderPath, 'dir') ~= 7
                return
            end

            items = dir(folderPath);
            for k = 1:numel(items)
                name = items(k).name;
                if ~items(k).isdir || strcmp(name, '.') || strcmp(name, '..')
                    continue
                end

                if ~isempty(regexp(name, '^\d{8}$', 'once')) && ~ismember(name, updateDates)
                    originalPath = fullfile(items(k).folder, name);
                    relativePath = fullfile(groupName, name);
                    hiddenPaths = app.hideSENUpdatePath(hiddenPaths, originalPath, hiddenRoot, relativePath);
                end
            end
        end

        function hiddenPaths = hideSENUpdatePath(app, hiddenPaths, originalPath, hiddenRoot, relativePath)
            if exist(originalPath, 'file') ~= 2 && exist(originalPath, 'dir') ~= 7
                return
            end

            hiddenPath = fullfile(hiddenRoot, relativePath);
            hiddenParent = fileparts(hiddenPath);
            if exist(hiddenParent, 'dir') ~= 7
                mkdir(hiddenParent);
            end

            movefile(originalPath, hiddenPath);
            hiddenPaths(end + 1).original = originalPath; %#ok<AGROW>
            hiddenPaths(end).hidden = hiddenPath;
        end

        function restoreSENUpdateWorkspace(app, hiddenPaths)
            for k = numel(hiddenPaths):-1:1
                hiddenPath = hiddenPaths(k).hidden;
                originalPath = hiddenPaths(k).original;

                if exist(hiddenPath, 'file') ~= 2 && exist(hiddenPath, 'dir') ~= 7
                    continue
                end

                originalParent = fileparts(originalPath);
                if exist(originalParent, 'dir') ~= 7
                    mkdir(originalParent);
                end

                if exist(originalPath, 'file') ~= 2 && exist(originalPath, 'dir') ~= 7
                    movefile(hiddenPath, originalPath);
                elseif exist(hiddenPath, 'dir') == 7 && exist(originalPath, 'dir') == 7
                    children = dir(hiddenPath);
                    for c = 1:numel(children)
                        childName = children(c).name;
                        if strcmp(childName, '.') || strcmp(childName, '..')
                            continue
                        end
                        childHidden = fullfile(children(c).folder, childName);
                        childOriginal = fullfile(originalPath, childName);
                        if exist(childOriginal, 'file') ~= 2 && exist(childOriginal, 'dir') ~= 7
                            movefile(childHidden, childOriginal);
                        end
                    end
                    if isempty(dir(fullfile(hiddenPath, '*')))
                        rmdir(hiddenPath);
                    end
                end
            end
            app.cleanupSENUpdateHiddenRoots(hiddenPaths);
        end

        function cleanupSENUpdateHiddenRoots(app, hiddenPaths)
            if isempty(hiddenPaths)
                return;
            end

            for k = 1:numel(hiddenPaths)
                hiddenPath = hiddenPaths(k).hidden;
                if exist(hiddenPath, 'file') == 2 || exist(hiddenPath, 'dir') == 7
                    return;
                end
            end

            hiddenRoots = {};
            for k = 1:numel(hiddenPaths)
                hiddenPath = hiddenPaths(k).hidden;
                parts = strsplit(hiddenPath, filesep);
                idx = find(strcmp(parts, '_update_hidden'), 1, 'first');
                if ~isempty(idx)
                    hiddenRoots{end + 1} = fullfile(parts{1:idx}); %#ok<AGROW>
                end
            end
            hiddenRoots = unique(hiddenRoots);

            for k = 1:numel(hiddenRoots)
                hiddenRoot = hiddenRoots{k};
                if exist(hiddenRoot, 'dir') == 7
                    try
                        rmdir(hiddenRoot, 's');
                    catch ME
                        updateOutput(app, sprintf('Warning: could not remove temporary update folder %s: %s', hiddenRoot, ME.message));
                    end
                end
            end
        end

        function hasDate = nameContainsAnySENDate(app, name, updateDates)
            hasDate = false;
            for k = 1:numel(updateDates)
                if contains(name, updateDates{k})
                    hasDate = true;
                    return
                end
            end
        end


        function run(app)


            if app.constellation == "SEN"

                condition = 0;

                while condition == 0

                    % Display a message in the UIAxes
                    updateOutput(app, 'Button pressed! Running the script...');

                    % Change lamp color to yellow (running)

                    try

                        %% SENTINEL-1 PREPROCESSING CODE ----------------------

                        pwd;
                        prep_folder = app.ProjectRoot;
                        addpath(genpath(app.InstallRoot));

                        %% ------------------ IMPORT OF THE REQUIRED VARIABLES --------------------


                        % DEFINE THE TYPE OF SLASH
                        par = filesep;

                        % READ OF THE INPUT VARIABLES FROM .MAT FILE
                        load(fullfile(phase_preprocessing_beta.dataFolder(app.ProjectRoot), 'input_preprocessing.mat'), 'python', ...
                                                'master_date', 'master_processing', 'auto_master', ...
                                                'polarisation', 'lon_min', 'lat_min', 'lon_max', 'lat_max', 'slaves_removal', 'dem_name', ...
                                                'dem_file', 'dem_name_coreg', 'dem_file_coreg', 'dem_resampling', 'first_step', 'coherence_tc', ...
                                                'epsg_code', 'gptbin_path', 'cpu', 'cache');

                        % PROJECT FOLDER
                        project_path_full = phase_preprocessing_beta.dataFolder(app.ProjectRoot);
                        if ~isfolder(fullfile(project_path_full,'snap2stamps','bin'))
                            mkdir(fullfile(project_path_full,'snap2stamps','bin'));
                        end

                        % GENERAL VARIABLES
                        space = (' ');
                        python = (python);
                        path_1_download = strcat(project_path_full, par, 'slaves');
                        update_processed_run = app.update_processed_data_SEN == 1;
                        update_new_dates = {};
                        update_hidden_paths = struct('original', {}, 'hidden', {});
                        update_cleanup = [];

                        % CONVERT VARIABLES FORMAT
                        if isa(master_date, 'datetime')
                            master_date = datestr(master_date, 'yyyymmdd');
                        end


                        %% ------------- CREATE FOLDERS ------------------

                        mkdir(project_path_full, 'slaves'); % create slaves folder
                        mkdir(project_path_full, 'master'); % create master folder

                        if update_processed_run
                            update_new_dates = app.getSENRootZipDates(project_path_full);
                            if isempty(update_new_dates)
                                error('Update already processed data is selected, but no new Sentinel-1 ZIP files were found in PHASE_Preprocessing/slaves.');
                            end

                            master_processing = 1;
                            first_step = 1;
                            slaves_removal = 1;
                            updateOutput(app, ['Update mode enabled. New date(s): ' strjoin(update_new_dates, ', ')]);
                            updateOutput(app, 'Update mode: master processing and slaves removal will be skipped.');
                            drawnow;
                        end

                        %% ------------------------ MASTER PROCESSING -----------------------------

                        if master_processing == 0 % check if the processing of the master has to be done or not

                            updateOutput(app, '----------------------- STEP 2: Master processing started -----------------------');

                            % 1. THE NUCLEAR RESET: Flatten all .zip files back to the root of /slaves

                            % Sweep 1: Rescue any Master .zip files from anywhere in /master
                            master_zips = dir(fullfile(project_path_full, 'master', '**', '*.zip'));
                            for i = 1:length(master_zips)
                                movefile(fullfile(master_zips(i).folder, master_zips(i).name), fullfile(project_path_full, 'slaves'));
                            end

                            % Sweep 2: Rescue all Slave .zip files from their date subfolders
                            slave_zips = dir(fullfile(project_path_full, 'slaves', '**', '*.zip'));
                            for i = 1:length(slave_zips)
                                % Only move it if it is NOT already in the root of /slaves
                                if ~strcmp(slave_zips(i).folder, fullfile(project_path_full, 'slaves'))
                                    movefile(fullfile(slave_zips(i).folder, slave_zips(i).name), fullfile(project_path_full, 'slaves'));
                                end
                            end

                            % Sweep 3: Delete all the now-empty date subfolders inside /slaves
                            slave_subdirs = dir(fullfile(project_path_full, 'slaves'));
                            for i = 1:length(slave_subdirs)
                                if slave_subdirs(i).isdir && ~strcmp(slave_subdirs(i).name, '.') && ~strcmp(slave_subdirs(i).name, '..')
                                    rmdir(fullfile(slave_subdirs(i).folder, slave_subdirs(i).name), 's');
                                end
                            end

                            % Clean old processing folders to prevent conflicts
                            folders_to_delete = {'split', 'coreg', 'ifg', 'lia', 'coherence', 'intensity'};
                            for i=1:length(folders_to_delete)
                                if exist(fullfile(project_path_full, folders_to_delete{i}), 'dir')
                                    rmdir(fullfile(project_path_full, folders_to_delete{i}), 's');
                                end
                                mkdir(fullfile(project_path_full, folders_to_delete{i}));
                            end

                            % Wipe master folder internals and recreate
                            if exist(fullfile(project_path_full, 'master'), 'dir')
                                rmdir(fullfile(project_path_full, 'master'), 's');
                            end
                            mkdir(fullfile(project_path_full, 'master'));

                            % 2. CREATE project_master.conf
                            projectfolder = strcat('PROJECTFOLDER=',project_path_full);
                            graphsfolder = strcat('GRAPHSFOLDER=', app.InstallRoot, par, 'PHASE_Preprocessing', par, 'snap2stamps', par, 'graphs');

                            if ispc
                                f_project_conf_master = fopen(strcat(project_path_full, '\snap2stamps\bin\project_master.conf'),'w');
                                j = 0;
                                while f_project_conf_master == -1 && j < 20
                                    f_project_conf_master = fopen(strcat(project_path_full, '\snap2stamps\bin\project_master.conf'),'w');
                                    j = j+1;
                                end
                                fprintf(f_project_conf_master,'%s\r\n', '######### CONFIGURATION FILE ######');
                                fprintf(f_project_conf_master,'%s\r\n', '# PROJECT DEFINITION');
                                fprintf(f_project_conf_master,'%s\r\n', projectfolder);
                                fprintf(f_project_conf_master,'%s\r\n', graphsfolder);
                                fprintf(f_project_conf_master,'%s\r\n', '# PROCESSING PARAMETERS');
                                fprintf(f_project_conf_master,'SWATHS=IW1,IW2,IW3\r\n');
                                fprintf(f_project_conf_master,'POLARISATION=%s\r\n', polarisation);
                                fprintf(f_project_conf_master,'AUTO_MASTER=%.0f\r\n', auto_master);
                                fprintf(f_project_conf_master,'MASTER_DATE=%s\r\n', master_date);
                                fprintf(f_project_conf_master,'%s\r\n', '# AOI BBOX DEFINITION');
                                fprintf(f_project_conf_master,'LONMIN=%.3f\r\n', lon_min);
                                fprintf(f_project_conf_master,'LATMIN=%.3f\r\n', lat_min);
                                fprintf(f_project_conf_master,'LONMAX=%.3f\r\n', lon_max);
                                fprintf(f_project_conf_master,'LATMAX=%.3f\r\n', lat_max);
                                fprintf(f_project_conf_master,'%s\r\n', '# SNAP GPT');
                                fprintf(f_project_conf_master,'GPTBIN_PATH=%s\r\n', gptbin_path);
                                fprintf(f_project_conf_master,'%s\r\n', '# COMPUTING RESOURCES');
                                fprintf(f_project_conf_master,'CPU=%.0f\r\n', cpu);
                                fprintf(f_project_conf_master,'CACHE=%s\r\n', cache);
                                fclose(f_project_conf_master);

                            elseif isunix && ~ismac
                                f_project_conf_master = fopen(strcat(project_path_full, '/snap2stamps/bin/project_master.conf'),'w');
                                j = 0;
                                while f_project_conf_master == -1 && j < 20
                                    f_project_conf_master = fopen(strcat(project_path_full, '/snap2stamps/bin/project_master.conf'),'w');
                                    j = j+1;
                                end
                                fprintf(f_project_conf_master,'%s\n', '######### CONFIGURATION FILE ######');
                                fprintf(f_project_conf_master,'%s\n', '# PROJECT DEFINITION');
                                fprintf(f_project_conf_master,'%s\n', projectfolder);
                                fprintf(f_project_conf_master,'%s\n', graphsfolder);
                                fprintf(f_project_conf_master,'%s\n', '# PROCESSING PARAMETERS');
                                fprintf(f_project_conf_master,'SWATHS=IW1,IW2,IW3\n');
                                fprintf(f_project_conf_master,'POLARISATION=%s\n', polarisation);
                                fprintf(f_project_conf_master,'AUTO_MASTER=%.0f\n', auto_master);
                                fprintf(f_project_conf_master,'MASTER_DATE=%s\n', master_date);
                                fprintf(f_project_conf_master,'%s\n', '# AOI BBOX DEFINITION');
                                fprintf(f_project_conf_master,'LONMIN=%.3f\n', lon_min);
                                fprintf(f_project_conf_master,'LATMIN=%.3f\n', lat_min);
                                fprintf(f_project_conf_master,'LONMAX=%.3f\n', lon_max);
                                fprintf(f_project_conf_master,'LATMAX=%.3f\n', lat_max);
                                fprintf(f_project_conf_master,'%s\n', '# SNAP GPT');
                                fprintf(f_project_conf_master,'GPTBIN_PATH=%s\n', gptbin_path);
                                fprintf(f_project_conf_master,'%s\n', '# COMPUTING RESOURCES');
                                fprintf(f_project_conf_master,'CPU=%.0f\n', cpu);
                                fprintf(f_project_conf_master,'CACHE=%s\n', cache);
                                fclose(f_project_conf_master);

                            elseif ismac
                                f_project_conf_master = fopen(strcat(project_path_full, '/snap2stamps/bin/project_master.conf'),'w');
                                j = 0;
                                while f_project_conf_master == -1 && j < 20
                                    f_project_conf_master = fopen(strcat(project_path_full, '/snap2stamps/bin/project_master.conf'),'w');
                                    j = j+1;
                                end
                                fprintf(f_project_conf_master,'%s\n', '######### CONFIGURATION FILE ######');
                                fprintf(f_project_conf_master,'%s\n', '# PROJECT DEFINITION');
                                fprintf(f_project_conf_master,'%s\n', projectfolder);
                                fprintf(f_project_conf_master,'%s\n', graphsfolder);
                                fprintf(f_project_conf_master,'%s\n', '# PROCESSING PARAMETERS');
                                fprintf(f_project_conf_master,'SWATHS=IW1,IW2,IW3\n');
                                fprintf(f_project_conf_master,'POLARISATION=%s\n', polarisation);
                                fprintf(f_project_conf_master,'AUTO_MASTER=%.0f\n', auto_master);
                                fprintf(f_project_conf_master,'MASTER_DATE=%s\n', master_date);
                                fprintf(f_project_conf_master,'%s\n', '# AOI BBOX DEFINITION');
                                fprintf(f_project_conf_master,'LONMIN=%.3f\n', lon_min);
                                fprintf(f_project_conf_master,'LATMIN=%.3f\n', lat_min);
                                fprintf(f_project_conf_master,'LONMAX=%.3f\n', lon_max);
                                fprintf(f_project_conf_master,'LATMAX=%.3f\n', lat_max);
                                fprintf(f_project_conf_master,'%s\n', '# SNAP GPT');
                                fprintf(f_project_conf_master,'GPTBIN_PATH=%s\n', gptbin_path);
                                fprintf(f_project_conf_master,'%s\n', '# COMPUTING RESOURCES');
                                fprintf(f_project_conf_master,'CPU=%.0f\n', cpu);
                                fprintf(f_project_conf_master,'CACHE=%s\n', cache);
                                fclose(f_project_conf_master);
                            end

                            % 3. EXECUTE MASTER SELECTOR AND MASTER PREP BATCH SCRIPT
                            updateOutput(app, 'Running automated Master Auto-Selector and Splitter...');
                            step_selector = phase_preprocessing_beta.scriptCommand(python, 'SEN_master_selector.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project_master.conf'), app.InstallRoot);
                            step_selector_cmd = step_selector;
                            step1_master = phase_preprocessing_beta.scriptCommand(python, 'SEN_splitting_master.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project_master.conf'), app.InstallRoot);
                            step_master_2 = step1_master;

                            path_1_master = fullfile(project_path_full, par, 'snap2stamps', par, 'bin');

                            if ispc
                                f_snap2stamps_master = fopen(strcat(project_path_full, '\snap2stamps\bin\snap2stamps_master.bat'),'w');
                                j = 0;
                                while f_snap2stamps_master == -1 && j < 20
                                    f_snap2stamps_master = fopen(strcat(project_path_full, '\snap2stamps\bin\snap2stamps_master.bat'),'w');
                                    j = j+1;
                                end
                                fprintf(f_snap2stamps_master,'@echo off \r\n');
                                fprintf(f_snap2stamps_master,'cd "%s"\r\n',path_1_master);
                                fprintf(f_snap2stamps_master,'%s\r\n',step_selector_cmd);
                                fprintf(f_snap2stamps_master,'%s\r\n',step_master_2);
                                fclose(f_snap2stamps_master);

                                path_2_master = fullfile(project_path_full, '\snap2stamps\bin\snap2stamps_master.bat');
                                phase_preprocessing_beta.runCommandLive(app, path_2_master, ...
                                    'Master selection and preparation', 3, 18, 1, 1);

                            elseif isunix && ~ismac
                                f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                j = 0;
                                while f_snap2stamps_master == -1 && j < 20
                                    f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                    j = j+1;
                                end
                                fprintf(f_snap2stamps_master,'#!/bin/bash \nset -e\n');
                                fprintf(f_snap2stamps_master,'cd "%s"\n',path_1_master);
                                fprintf(f_snap2stamps_master,'%s\n',step_selector_cmd);
                                fprintf(f_snap2stamps_master,'%s\n',step_master_2);
                                fclose(f_snap2stamps_master);

                                path_2_master = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'));
                                chmod = ['chmod +x' space path_2_master];
                                phase_preprocessing_beta.runCommandLive(app, path_2_master, ...
                                    'Master selection and preparation', 3, 18, 1, 1);

                            elseif ismac
                                f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                j = 0;
                                while f_snap2stamps_master == -1 && j < 20
                                    f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                    j = j+1;
                                end
                                fprintf(f_snap2stamps_master,'#!/bin/bash \nset -e\n');
                                fprintf(f_snap2stamps_master,'cd "%s"\n',path_1_master);
                                fprintf(f_snap2stamps_master,'%s\n',step_selector_cmd);
                                fprintf(f_snap2stamps_master,'%s\n',step_master_2);
                                fclose(f_snap2stamps_master);

                                path_2_master = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'));
                                chmod = ['chmod +x' space path_2_master];
                                phase_preprocessing_beta.runCommandLive(app, path_2_master, ...
                                    'Master selection and preparation', 3, 18, 1, 1);
                            end

                            % FLATTEN DIRECTORY: Move files up one level to clean the structure
                            master_subfolders = dir(fullfile(project_path_full, 'master', 'S1*'));
                            for idx = 1:length(master_subfolders)
                                if master_subfolders(idx).isdir
                                    subfolder_path = fullfile(master_subfolders(idx).folder, master_subfolders(idx).name);

                                    % Move all contents (.zip, .dim, .data) up to /master
                                    all_contents = dir(fullfile(subfolder_path, '*'));
                                    for j = 1:length(all_contents)
                                        if ~strcmp(all_contents(j).name, '.') && ~strcmp(all_contents(j).name, '..')
                                            movefile(fullfile(all_contents(j).folder, all_contents(j).name), fullfile(project_path_full, 'master'));
                                        end
                                    end
                                    % Delete the now-empty 67-character folder
                                    rmdir(subfolder_path, 's');
                                end
                            end

                            % 4. VALIDATE OUTPUT
                            master_dim_struct = dir(fullfile(project_path_full, 'master', '*_split_*_Orb.dim'));

                            if isempty(master_dim_struct)
                                msg = 'ERROR: the master image has not been correctly processed. Check the snap2stamps logs.';
                                error(msg);
                            end
                            master_file_destination = fullfile(master_dim_struct(1).folder, master_dim_struct(1).name);

                            updateOutput(app, '----------------------- STEP 2: Master processing finished -----------------------');

                        else
                            % IF MASTER PROCESSING IS SKIPPED, DYNAMICALLY FIND THE EXISTING MASTER
                            master_dim_struct = dir(fullfile(project_path_full, 'master', '*_split_*_Orb.dim'));

                            if isempty(master_dim_struct)
                                msg = 'ERROR: No processed master found in /master folder! You must run Master Processing first.';
                                error(msg);
                            end
                            master_file_destination = fullfile(master_dim_struct(1).folder, master_dim_struct(1).name);

                            updateOutput(app, '----------------------- STEP 2: Master processing skipped -----------------------');
                            drawnow;
                        end

                        if app.StopFlag
                            break;
                        end


                        %% ------------------------ SLAVES PROCESSING -----------------------------

                        if update_processed_run
                            updateOutput(app, 'Update mode: preparing temporary workspace for new slave images.');
                            drawnow;
                            update_hidden_paths = app.isolateSENUpdateWorkspace(project_path_full, update_new_dates);
                            update_cleanup = onCleanup(@() app.restoreSENUpdateWorkspace(update_hidden_paths));
                            updateOutput(app, 'Update mode: existing slave products temporarily hidden from the slaves pipeline.');
                            drawnow;
                        end

                        % CREATE THE PROJECT.CONF FILE FOR SNAP2STAMPS SLAVES PROCESSING
                        projectfolder = strcat('PROJECTFOLDER=',project_path_full); % path of the project folder
                        graphsfolder = strcat('GRAPHSFOLDER=', app.InstallRoot, par, 'PHASE_Preprocessing', par, 'snap2stamps', par, 'graphs'); % path of the graphs folder of snap2stamps
                        masterfolder = strcat('MASTER=', master_file_destination); % master image path

                        if ispc
                            f_project_conf_slaves = fopen(strcat(project_path_full, '\snap2stamps\bin\project.conf'),'w');
                            j = 0;
                            while f_project_conf_slaves == -1 && j < 20
                                f_project_conf_slaves = fopen(strcat(project_path_full, '\snap2stamps\bin\project.conf'),'w');
                                j = j+1;
                            end
                            fprintf(f_project_conf_slaves,'%s\r\n', '######### CONFIGURATION FILE ######');
                            fprintf(f_project_conf_slaves,'%s\r\n', '###################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# PROJECT DEFINITION');
                            fprintf(f_project_conf_slaves,'%s\r\n', projectfolder);
                            fprintf(f_project_conf_slaves,'%s\r\n', graphsfolder);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# PROCESSING PARAMETERS');
                            fprintf(f_project_conf_slaves,'SWATHS=IW1,IW2,IW3\r\n');
                            fprintf(f_project_conf_slaves,'%s\r\n', masterfolder);
                            fprintf(f_project_conf_slaves,'POLARISATION=%s\r\n', polarisation);
                            fprintf(f_project_conf_slaves,'TC_COHERENCE=%.0f\n', coherence_tc);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# DEM DEFINITION');
                            fprintf(f_project_conf_slaves,'DEMNAME=%s\r\n', dem_name);
                            fprintf(f_project_conf_slaves,'DEMFILE=%s\r\n', dem_file);
                            fprintf(f_project_conf_slaves,'COREGDTMNAME=%s\r\n', dem_name_coreg);
                            fprintf(f_project_conf_slaves,'COREGDTMFILE=%s\r\n', dem_file_coreg);
                            fprintf(f_project_conf_slaves,'DEMRESAMPLING=%s\r\n', dem_resampling);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# AOI BBOX DEFINITION');
                            fprintf(f_project_conf_slaves,'LONMIN=%.3f\r\n', lon_min);
                            fprintf(f_project_conf_slaves,'LATMIN=%.3f\r\n', lat_min);
                            fprintf(f_project_conf_slaves,'LONMAX=%.3f\r\n', lon_max);
                            fprintf(f_project_conf_slaves,'LATMAX=%.3f\r\n', lat_max);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# SNAP GPT');
                            fprintf(f_project_conf_slaves,'GPTBIN_PATH=%s\r\n', gptbin_path);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# COMPUTING RESOURCES TO EMPLOY');
                            fprintf(f_project_conf_slaves,'CPU=%.0f\r\n', cpu);
                            fprintf(f_project_conf_slaves,'CACHE=%s\r\n', cache);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            f_project_conf_slaves = fclose(f_project_conf_slaves);
                        elseif isunix && ~ismac
                            f_project_conf_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/project.conf'),'w');
                            j = 0;
                            while f_project_conf_slaves == -1 && j < 20
                                f_project_conf_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/project.conf'),'w');
                                j = j+1;
                            end
                            fprintf(f_project_conf_slaves,'%s\n', '######### CONFIGURATION FILE ######');
                            fprintf(f_project_conf_slaves,'%s\n', '###################################');
                            fprintf(f_project_conf_slaves,'%s\n', '# PROJECT DEFINITION');
                            fprintf(f_project_conf_slaves,'%s\n', projectfolder);
                            fprintf(f_project_conf_slaves,'%s\n', graphsfolder);
                            fprintf(f_project_conf_slaves,'%s\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\n', '# PROCESSING PARAMETERS');
                            fprintf(f_project_conf_slaves,'SWATHS=IW1,IW2,IW3\r\n');
                            fprintf(f_project_conf_slaves,'%s\r\n', masterfolder);
                            fprintf(f_project_conf_slaves,'POLARISATION=%s\r\n', polarisation);
                            fprintf(f_project_conf_slaves,'TC_COHERENCE=%.0f\n', coherence_tc);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# DEM DEFINITION');
                            fprintf(f_project_conf_slaves,'DEMNAME=%s\r\n', dem_name);
                            fprintf(f_project_conf_slaves,'DEMFILE=%s\r\n', dem_file);
                            fprintf(f_project_conf_slaves,'COREGDTMNAME=%s\r\n', dem_name_coreg);
                            fprintf(f_project_conf_slaves,'COREGDTMFILE=%s\r\n', dem_file_coreg);
                            fprintf(f_project_conf_slaves,'DEMRESAMPLING=%s\r\n', dem_resampling);
                            fprintf(f_project_conf_slaves,'%s\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\n', '# AOI BBOX DEFINITION');
                            fprintf(f_project_conf_slaves,'LONMIN=%.3f\n', lon_min);
                            fprintf(f_project_conf_slaves,'LATMIN=%.3f\n', lat_min);
                            fprintf(f_project_conf_slaves,'LONMAX=%.3f\n', lon_max);
                            fprintf(f_project_conf_slaves,'LATMAX=%.3f\n', lat_max);
                            fprintf(f_project_conf_slaves,'%s\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\n', '# SNAP GPT');
                            fprintf(f_project_conf_slaves,'GPTBIN_PATH=%s\n', gptbin_path);
                            fprintf(f_project_conf_slaves,'%s\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\n', '# COMPUTING RESOURCES TO EMPLOY');
                            fprintf(f_project_conf_slaves,'CPU=%.0f\n', cpu);
                            fprintf(f_project_conf_slaves,'CACHE=%s\n', cache);
                            fprintf(f_project_conf_slaves,'%s\n', '##################################');
                            f_project_conf_slaves = fclose(f_project_conf_slaves);
                        elseif ismac
                            f_project_conf_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/project.conf'),'w');
                            j = 0;
                            while f_project_conf_slaves == -1 && j < 20
                                f_project_conf_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/project.conf'),'w');
                                j = j+1;
                            end
                            fprintf(f_project_conf_slaves,'%s\n', '######### CONFIGURATION FILE ######');
                            fprintf(f_project_conf_slaves,'%s\n', '###################################');
                            fprintf(f_project_conf_slaves,'%s\n', '# PROJECT DEFINITION');
                            fprintf(f_project_conf_slaves,'%s\n', projectfolder);
                            fprintf(f_project_conf_slaves,'%s\n', graphsfolder);
                            fprintf(f_project_conf_slaves,'%s\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\n', '# PROCESSING PARAMETERS');
                            fprintf(f_project_conf_slaves,'SWATHS=IW1,IW2,IW3\r\n');
                            fprintf(f_project_conf_slaves,'%s\r\n', masterfolder);
                            fprintf(f_project_conf_slaves,'POLARISATION=%s\r\n', polarisation);
                            fprintf(f_project_conf_slaves,'TC_COHERENCE=%.0f\n', coherence_tc);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# DEM DEFINITION');
                            fprintf(f_project_conf_slaves,'DEMNAME=%s\r\n', dem_name);
                            fprintf(f_project_conf_slaves,'DEMFILE=%s\r\n', dem_file);
                            fprintf(f_project_conf_slaves,'COREGDTMNAME=%s\r\n', dem_name_coreg);
                            fprintf(f_project_conf_slaves,'COREGDTMFILE=%s\r\n', dem_file_coreg);
                            fprintf(f_project_conf_slaves,'DEMRESAMPLING=%s\r\n', dem_resampling);
                            fprintf(f_project_conf_slaves,'%s\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\n', '# AOI BBOX DEFINITION');
                            fprintf(f_project_conf_slaves,'LONMIN=%.3f\n', lon_min);
                            fprintf(f_project_conf_slaves,'LATMIN=%.3f\n', lat_min);
                            fprintf(f_project_conf_slaves,'LONMAX=%.3f\n', lon_max);
                            fprintf(f_project_conf_slaves,'LATMAX=%.3f\n', lat_max);
                            fprintf(f_project_conf_slaves,'%s\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\n', '# SNAP GPT');
                            fprintf(f_project_conf_slaves,'GPTBIN_PATH=%s\n', gptbin_path);
                            fprintf(f_project_conf_slaves,'%s\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\n', '# COMPUTING RESOURCES TO EMPLOY');
                            fprintf(f_project_conf_slaves,'CPU=%.0f\n', cpu);
                            fprintf(f_project_conf_slaves,'CACHE=%s\n', cache);
                            fprintf(f_project_conf_slaves,'%s\n', '##################################');
                            f_project_conf_slaves = fclose(f_project_conf_slaves);
                        end

                        % EXECUTE THE .CONF FILE VIA A BATCH/BASH FILE
                        dp = phase_preprocessing_beta.skipPrefix();
                        step1_slaves = phase_preprocessing_beta.scriptCommand(python, 'SEN_slaves_prep.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % slaves preparation
                        step2_slaves = phase_preprocessing_beta.scriptCommand(python, 'SEN_splitting_slaves.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % slaves splitting & apply orbits
                        step3_slaves = phase_preprocessing_beta.scriptCommand(python, 'SEN_coreg_ifg_topsar.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % coregistration & interferogram
                        step4_slaves = phase_preprocessing_beta.scriptCommand(python, 'SEN_stamps_export.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % StaMPS export
                        step5_slaves = phase_preprocessing_beta.scriptCommand(python, 'SEN_average_intensity.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % average instensity
                        step6_slaves = phase_preprocessing_beta.scriptCommand(python, 'SEN_terrain_correction.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % terrain corrected coherence and lia

                            % CASES DEFENDING ON FIRST STEP
                            if isnumeric(first_step)
                                first_step_num = first_step;
                            elseif ischar(first_step)
                                first_step_num = str2double(first_step);
                            end
                            if first_step_num == 1
                                step_slaves_1 = step1_slaves;
                                step_slaves_2 = step2_slaves;
                                step_slaves_3 = step3_slaves;
                                step_slaves_4 = step4_slaves;
                                step_slaves_5 = step5_slaves;
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                end
                            elseif first_step_num == 2
                                step_slaves_1 = [dp step1_slaves];
                                step_slaves_2 = step2_slaves;
                                step_slaves_3 = step3_slaves;
                                step_slaves_4 = step4_slaves;
                                step_slaves_5 = step5_slaves;
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                end
                            elseif first_step_num == 3
                                step_slaves_1 = [dp step1_slaves];
                                step_slaves_2 = [dp step2_slaves];
                                step_slaves_3 = step3_slaves;
                                step_slaves_4 = step4_slaves;
                                step_slaves_5 = step5_slaves;
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                end
                            elseif first_step_num == 4
                                step_slaves_1 = [dp step1_slaves];
                                step_slaves_2 = [dp step2_slaves];
                                step_slaves_3 = [dp step3_slaves];
                                step_slaves_4 = step4_slaves;
                                step_slaves_5 = step5_slaves;
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                end
                             elseif first_step_num == 5
                                step_slaves_1 = [dp step1_slaves];
                                step_slaves_2 = [dp step2_slaves];
                                step_slaves_3 = [dp step3_slaves];
                                step_slaves_4 = [dp step4_slaves];
                                step_slaves_5 = step5_slaves;
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                end
                            elseif first_step_num == 6
                                step_slaves_1 = [dp step1_slaves];
                                step_slaves_2 = [dp step2_slaves];
                                step_slaves_3 = [dp step3_slaves];
                                step_slaves_4 = [dp step4_slaves];
                                step_slaves_5 = [dp step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                else
                                    updateOutput(app, ['You are running just the terrain correction for coherence and LIA bands, but you have set tc_coherence = 1. ' ...
                                        'Please change it to 0 to perform this step.'])
                                end
                            end

                        path_1_slaves = fullfile(project_path_full, par, 'snap2stamps', par, 'bin');
                        slaves_proof = strcat(project_path_full, par, 'slaves', par, 'dummyDownload.txt');

                        if ispc
                            f_snap2stamps_slaves = fopen(strcat(project_path_full, '\snap2stamps\bin\snap2stamps_slaves.bat'),'w');
                            j = 0;
                            while f_snap2stamps_slaves == -1 && j < 20
                                f_snap2stamps_slaves = fopen(strcat(project_path_full, '\snap2stamps\bin\snap2stamps_slaves.bat'),'w');
                                j = j+1;
                            end
                            fprintf(f_snap2stamps_slaves,'@echo off \r\n');
                            fprintf(f_snap2stamps_slaves,'cd "%s"\r\n',path_1_slaves);
                            fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_1);
                            fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_2);
                            fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_3);
                            fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_4);
                            fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_5);
                            if coherence_tc == 0
                                fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_6);
                            end
                            fprintf(f_snap2stamps_slaves,'cd "%s"\r\n',path_1_download);
                            fprintf(f_snap2stamps_slaves,'type nul > dummyDownload.txt \r\n');
                            fprintf(f_snap2stamps_slaves,'exit \r\n');
                            f_snap2stamps_slaves = fclose(f_snap2stamps_slaves);
                        elseif isunix && ~ismac
                            f_snap2stamps_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh'),'w');
                            j = 0;
                            while f_snap2stamps_slaves == -1 && j < 20
                                f_snap2stamps_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh'),'w');
                                j = j+1;
                            end
                            fprintf(f_snap2stamps_slaves,'#!/bin/bash \nset -e\n');
                            fprintf(f_snap2stamps_slaves,'cd "%s"\n',path_1_slaves);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_1);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_2);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_3);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_4);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_5);
                            if coherence_tc == 0
                                fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_6);
                            end
                            fprintf(f_snap2stamps_slaves,'sleep 5 \n');
                            fprintf(f_snap2stamps_slaves,'cd "%s"\n',path_1_download);
                            fprintf(f_snap2stamps_slaves,'touch dummyDownload.txt \n');
                            f_snap2stamps_slaves = fclose(f_snap2stamps_slaves);
                        elseif ismac
                            f_snap2stamps_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh'),'w');
                            j = 0;
                            while f_snap2stamps_slaves == -1 && j < 20
                                f_snap2stamps_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh'),'w');
                                j = j+1;
                            end
                            fprintf(f_snap2stamps_slaves,'#!/bin/bash \nset -e\n');
                            fprintf(f_snap2stamps_slaves,'cd "%s"\n',path_1_slaves);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_1);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_2);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_3);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_4);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_5);
                            if coherence_tc == 0
                                fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_6);
                            end
                            fprintf(f_snap2stamps_slaves,'sleep 5 \n');
                            fprintf(f_snap2stamps_slaves,'cd "%s"\n',path_1_download);
                            fprintf(f_snap2stamps_slaves,'touch dummyDownload.txt \n');
                            f_snap2stamps_slaves = fclose(f_snap2stamps_slaves);
                        end

                        updateOutput(app, '----------------------- STEP 3: Slaves processing started -----------------------');
                        drawnow;

                        if ispc
                            path_2_slaves = fullfile(project_path_full, '\snap2stamps\bin\snap2stamps_slaves.bat');
                            phase_preprocessing_beta.runCommandLive(app, path_2_slaves, ...
                                'Slave processing pipeline', 18, 92, first_step_num, 6);
                        elseif isunix && ~ismac
                            path_2_slaves = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh'));
                            chmod_in = ('chmod +x');
                            chmod = [chmod_in space path_2_slaves];

                            phase_preprocessing_beta.runCommandLive(app, path_2_slaves, ...
                                'Slave processing pipeline', 18, 92, first_step_num, 6);
                        elseif ismac
                            path_2_slaves = strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh');
                            chmod_in = 'chmod +x';
                            chmod = [chmod_in space path_2_slaves];
                            % Executed without opening Terminal by runCommandLive.

                            phase_preprocessing_beta.runCommandLive(app, path_2_slaves, ...
                                'Slave processing pipeline', 18, 92, first_step_num, 6);
                        end

                        % wait until the cmd has finished the work
                        while exist(slaves_proof,'file') == 0
                              pause(1);
                        end
                        delete(slaves_proof);

                        if update_processed_run
                            app.restoreSENUpdateWorkspace(update_hidden_paths);
                            update_hidden_paths = struct('original', {}, 'hidden', {});
                            clear update_cleanup
                            updateOutput(app, 'Update mode: existing slave products restored after slaves pipeline.');

                            average_intensity_proof = strcat(project_path_full, par, 'slaves', par, 'dummyAverageIntensity.txt');
                            app.deleteFileIfExists(average_intensity_proof);
                            updateOutput(app, 'Update mode: full-stack average intensity started.');

                            if ispc
                                f_average_intensity = fopen(strcat(project_path_full, '\snap2stamps\bin\snap2stamps_update_average_intensity.bat'),'w');
                                j = 0;
                                while f_average_intensity == -1 && j < 20
                                    f_average_intensity = fopen(strcat(project_path_full, '\snap2stamps\bin\snap2stamps_update_average_intensity.bat'),'w');
                                    j = j+1;
                                end
                                fprintf(f_average_intensity,'@echo off \r\n');
                                fprintf(f_average_intensity,'cd "%s"\r\n',path_1_slaves);
                                fprintf(f_average_intensity,'%s\r\n',step5_slaves);
                                fprintf(f_average_intensity,'cd "%s"\r\n',path_1_download);
                                fprintf(f_average_intensity,'type nul > dummyAverageIntensity.txt \r\n');
                                fprintf(f_average_intensity,'exit \r\n');
                                f_average_intensity = fclose(f_average_intensity);
                                path_2_average_intensity = fullfile(project_path_full, '\snap2stamps\bin\snap2stamps_update_average_intensity.bat');
                                phase_preprocessing_beta.runCommandLive(app, path_2_average_intensity, ...
                                    'Full-stack average intensity', 78, 92, 5, 5);
                            elseif isunix && ~ismac
                                f_average_intensity = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_update_average_intensity.sh'),'w');
                                j = 0;
                                while f_average_intensity == -1 && j < 20
                                    f_average_intensity = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_update_average_intensity.sh'),'w');
                                    j = j+1;
                                end
                                fprintf(f_average_intensity,'#!/bin/bash \nset -e\n');
                                fprintf(f_average_intensity,'cd "%s"\n',path_1_slaves);
                                fprintf(f_average_intensity,'%s\n',step5_slaves);
                                fprintf(f_average_intensity,'sleep 5 \n');
                                fprintf(f_average_intensity,'cd "%s"\n',path_1_download);
                                fprintf(f_average_intensity,'touch dummyAverageIntensity.txt \n');
                                f_average_intensity = fclose(f_average_intensity);
                                path_2_average_intensity = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_update_average_intensity.sh'));
                                chmod_in = ('chmod +x');
                                chmod = [chmod_in space path_2_average_intensity];
                                xterm = ('sudo xterm -e');

                                phase_preprocessing_beta.runCommandLive(app, path_2_average_intensity, ...
                                    'Full-stack average intensity', 78, 92, 5, 5);
                            elseif ismac
                                f_average_intensity = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_update_average_intensity.sh'),'w');
                                j = 0;
                                while f_average_intensity == -1 && j < 20
                                    f_average_intensity = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_update_average_intensity.sh'),'w');
                                    j = j+1;
                                end
                                fprintf(f_average_intensity,'#!/bin/bash \nset -e\n');
                                fprintf(f_average_intensity,'cd "%s"\n',path_1_slaves);
                                fprintf(f_average_intensity,'%s\n',step5_slaves);
                                fprintf(f_average_intensity,'sleep 5 \n');
                                fprintf(f_average_intensity,'cd "%s"\n',path_1_download);
                                fprintf(f_average_intensity,'touch dummyAverageIntensity.txt \n');
                                f_average_intensity = fclose(f_average_intensity);
                                path_2_average_intensity = strcat(project_path_full, '/snap2stamps/bin/snap2stamps_update_average_intensity.sh');
                                chmod_in = 'chmod +x';
                                chmod = [chmod_in space path_2_average_intensity];
                                % Executed without opening Terminal by runCommandLive.

                                phase_preprocessing_beta.runCommandLive(app, path_2_average_intensity, ...
                                    'Full-stack average intensity', 78, 92, 5, 5);
                            end

                            while exist(average_intensity_proof,'file') == 0
                                pause(1);
                            end
                            delete(average_intensity_proof);
                            updateOutput(app, 'Update mode: full-stack average intensity finished.');
                        end

                        % DYNAMIC WAIT FOR GEO-TIFF CREATION (Timeout: 5 minutes)
                        %updateOutput(app, 'Waiting for .tif files to be written on disk...');
                        timeout = 300; % Maximum wait time in seconds
                        elapsed = 0;
                        intensity_folder = fullfile(project_path_full, 'intensity');

                        % Wait until the folder exists and is not empty, but do not exceed the timeout
                        while (~exist(intensity_folder, 'dir') || isempty(dir(fullfile(intensity_folder, '*.tif')))) && elapsed < timeout
                            pause(5); % Check every 5 seconds
                            elapsed = elapsed + 5;
                        end

                        if elapsed >= timeout
                            updateOutput(app, 'WARNING: Timeout exceeded. The .tif files were not found.');
                        else
                            %updateOutput(app, '.tif files found, starting reprojection...');
                            % Leave a 2-second safety margin for the OS to release the file lock
                            pause(2);
                        end

                        % APPLY THE CORRECT PROJECTION TO THE AVERAGE INTENSITY GEOTIFF
                        if first_step_num <= 5
                            intensity_folder = fullfile(project_path_full, 'intensity');
                            if exist(intensity_folder, 'dir')
                                intensity_files = dir(fullfile(intensity_folder, '*.tif'));
                                for i = 1:length(intensity_files)
                                    file_path = fullfile(intensity_folder, intensity_files(i).name);
                                    app.reprojectGeoTiffIfGeographic(file_path, epsg_code);
                                end
                            end
                        end

                        coherence_folder = fullfile(project_path_full, 'coherence');
                        if exist(coherence_folder, 'dir')
                            coherence_files = dir(fullfile(coherence_folder, '*.tif'));
                            for i = 1:length(coherence_files)
                                file_path = fullfile(coherence_folder, coherence_files(i).name);
                                app.reprojectGeoTiffIfGeographic(file_path, epsg_code);
                            end
                        end

                        lia_folder = fullfile(project_path_full, 'lia');
                        if exist(lia_folder, 'dir')
                            lia_files = dir(fullfile(lia_folder, '*.tif'));
                            for i = 1:length(lia_files)
                                file_path = fullfile(lia_folder, lia_files(i).name);
                                app.reprojectGeoTiffIfGeographic(file_path, epsg_code, 2);
                            end
                        end

                        updateOutput(app, '----------------------- STEP 3: Slaves processing finished -----------------------');

                        if app.StopFlag
                            break;
                        end

                        % DYNAMIC STAMPS FOLDER AUTOMATION (ORBIT & DATE EXTRACTION)
                        coreg_dir = fullfile(project_path_full, 'coreg');
                        dim_files = dir(fullfile(coreg_dir, '*.dim'));

                        if isempty(dim_files)
                            updateOutput(app, 'ERROR: No coregistered files found to extract dates/orbit!');
                            return;
                        end

                        % 1. Extract all dates from the filenames to find absolute Min and Max
                        all_dates = NaT(0); % Initialize empty datetime array
                        for idx = 1:length(dim_files)
                            % Parse the YYYYMMDD_YYYYMMDD.dim filename structure
                            tokens = regexp(dim_files(idx).name, '(\d{8})_(\d{8})\.dim', 'tokens');
                            if ~isempty(tokens)
                                all_dates(end+1, 1) = datetime(tokens{1}{1}, 'InputFormat', 'yyyyMMdd');
                                all_dates(end+1, 1) = datetime(tokens{1}{2}, 'InputFormat', 'yyyyMMdd');
                            end
                        end

                        % Sort and format to English MmmYY (e.g., Jul19, Aug23)
                        min_date = min(all_dates);
                        max_date = max(all_dates);
                        date_str_start = char(datetime(min_date, 'Format', 'MMMyy', 'Locale', 'en_US'));
                        date_str_end   = char(datetime(max_date, 'Format', 'MMMyy', 'Locale', 'en_US'));

                        % 2. Read the Ground Truth Orbit from the Master .dim XML
                        dim_path = fullfile(dim_files(1).folder, dim_files(1).name);
                        dim_text = fileread(dim_path);

                        orbit_type = 'UNKNOWN';
                        if contains(upper(dim_text), 'ASCENDING')
                            orbit_type = 'ASC';
                        elseif contains(upper(dim_text), 'DESCENDING')
                            orbit_type = 'DES';
                        end

                        % 3. Automatically Construct the Target Folder Name!
                        stamps_folder = sprintf('%s_%s_%s', orbit_type, date_str_start, date_str_end);
                        updateOutput(app, ['Dynamically generated StaMPS folder name: ' stamps_folder]);

                        project_parent_path_full = phase_preprocessing_beta.stampsFolder(app.ProjectRoot);
                        stamps_folder_full = fullfile(project_parent_path_full, stamps_folder);
                        if ~isfolder(stamps_folder_full)
                            mkdir(stamps_folder_full);
                            updateOutput(app, ['Created StaMPS dataset folder: ' stamps_folder_full]);
                        end
                        % The beta launches the editable StaMPS module from its
                        % canonical installation. No PHASE_StaMPS.mlapp is copied
                        % or required at runtime.
                        stamps_diagnostic = fullfile(app.InstallRoot, 'PHASE_Preprocessing', 'diagnose_PHASE_StaMPS.m');
                        if isfile(stamps_diagnostic)
                            copyfile(stamps_diagnostic, stamps_folder_full, 'f');
                        end

                        % EVENTUAL REMOVAL OF THE SLAVES IMAGES TO SAVE DISK SPACE
                        if slaves_removal == 0 % check if the removal of the slaves images has to be done or not

                            updateOutput(app, '----------------------- STEP 4: Slaves removal started -----------------------');

                            rmdir(strcat(project_path_full, par, 'slaves'), 's'); % removal of slaves folder to save disk space
                            mkdir(project_path_full, 'slaves'); % create slaves folder

                            updateOutput(app, '----------------------- STEP 4: Slaves removal finished -----------------------');

                            else

                            updateOutput(app, '----------------------- STEP 4: Slaves removal skipped -----------------------');

                        end

                        % Seed a new dataset with the installed configuration, but never
                        % overwrite an existing input_StaMPS.mat: it may contain parameters
                        % the user tuned specifically for this processing run.
                        src_input_mat = fullfile(project_path_full, 'input_StaMPS.mat');
                        dst_input_mat = fullfile(stamps_folder_full, 'input_StaMPS.mat');
                        if isfile(dst_input_mat)
                            updateOutput(app, 'Preserved the existing dataset input_StaMPS.mat.');
                        elseif isfile(src_input_mat)
                            copyfile(src_input_mat, dst_input_mat);
                            updateOutput(app, 'Copied input_StaMPS.mat to the new StaMPS processing folder.');
                        else
                            updateOutput(app, ['NOTICE: input_StaMPS.mat is not available yet. ', ...
                                'PHASE_StaMPS_beta will create it when the initial settings are saved.']);
                        end

                        stamps_app_full = fullfile(project_parent_path_full, stamps_folder);
                        updateOutput(app, ['Preprocessing completed. StaMPS dataset folder: ' stamps_app_full]);
                        updateOutput(app, 'Open the StaMPS PSI section in the PHASE hub when ready.');

                        %% ----------------------------------------------------

                        % Display a message when the script is done
                        updateOutput(app, 'Script SEN_Preprocessing.m completed.');

                        % Change lamp color to green (success)

                    catch ME

                        % If an error occurs
                        if strcmp(ME.identifier,'PHASE:ProcessingStopped')
                            updateOutput(app, 'Preprocessing was force-stopped by the user.');
                        else
                            updateOutput(app, 'An error occurred during script execution. Please check each step log file!');
                        end

                        % Change lamp color to red (error)

                        % Show the actual error
                        rethrow(ME)

                    end

                    condition = 1;
                end

            elseif app.constellation == "CSK"

                condition = 0;

                while condition == 0


                    % Display a message in the UIAxes
                    updateOutput(app, 'Button pressed! Running the script...');
                    drawnow;

                    % Change lamp color to yellow (running)

                    try

                        %% COSMO-SKYMED PROCESSING ----------------------------

                        pwd;
                        prep_folder = app.ProjectRoot;
                        addpath(genpath(app.InstallRoot));

                        %% ------------------ IMPORT OF THE REQUIRED VARIABLES --------------------

                        % READ OF THE INPUT VARIABLES
                        par = filesep;

                        load(fullfile(phase_preprocessing_beta.dataFolder(app.ProjectRoot), 'input_preprocessing.mat'), 'python', ...
                                                'master_date', 'auto_master', 'master_processing', 'lon_min', ...
                                                'lat_min', 'lon_max', 'lat_max', 'slaves_removal', 'dem_name', 'dem_file', ...
                                                'first_step', 'num_gcp', 'coherence_tc', 'epsg_code', 'gptbin_path', 'cpu', 'cache');

                        % PROJECT FOLDER
                        project_path_full = phase_preprocessing_beta.dataFolder(app.ProjectRoot);
                        if ~isfolder(fullfile(project_path_full,'snap2stamps','bin'))
                            mkdir(fullfile(project_path_full,'snap2stamps','bin'));
                        end
                        mkdir(project_path_full, 'slaves'); % create slaves folder
                        mkdir(project_path_full, 'master'); % create master folder

                        % GENERAL VARIABLES
                        space = (' ');
                        python = (python);
                        path_1_download = strcat(project_path_full, par, 'slaves');

                        % CONVERT VARIABLES FORMAT
                        if isa(master_date, 'datetime')
                            master_date = datestr(master_date, 'yyyymmdd');
                        end


                        %% ------------------------ MASTER PROCESSING -----------------------------

                        if master_processing == 0 % check if the processing of the master has to be done or not

                            updateOutput(app, '----------------------- STEP 1: Master processing started -----------------------');

                            % 1. THE NUCLEAR RESET: Flatten all .h5 files back to the root of /slaves

                            % Sweep 1: Rescue any Master .h5 files from anywhere in /master
                            master_h5 = dir(fullfile(project_path_full, 'master', '**', '*.h5'));
                            for i = 1:length(master_h5)
                                movefile(fullfile(master_h5(i).folder, master_h5(i).name), fullfile(project_path_full, 'slaves'));
                            end

                            % Sweep 2: Rescue all Slave .h5 files from their date subfolders
                            slave_h5 = dir(fullfile(project_path_full, 'slaves', '**', '*.h5'));
                            for i = 1:length(slave_h5)
                                % Only move it if it is NOT already in the root of /slaves
                                if ~strcmp(slave_h5(i).folder, fullfile(project_path_full, 'slaves'))
                                    movefile(fullfile(slave_h5(i).folder, slave_h5(i).name), fullfile(project_path_full, 'slaves'));
                                end
                            end

                            % Sweep 3: Delete all the now-empty date subfolders inside /slaves
                            slave_subdirs = dir(fullfile(project_path_full, 'slaves'));
                            for i = 1:length(slave_subdirs)
                                if slave_subdirs(i).isdir && ~strcmp(slave_subdirs(i).name, '.') && ~strcmp(slave_subdirs(i).name, '..')
                                    rmdir(fullfile(slave_subdirs(i).folder, slave_subdirs(i).name), 's');
                                end
                            end

                            % Clean old processing folders to prevent conflicts
                            folders_to_delete = {'subset', 'coreg', 'ifg', 'lia', 'coherence', 'intensity'};
                            for i=1:length(folders_to_delete)
                                if exist(fullfile(project_path_full, folders_to_delete{i}), 'dir')
                                    rmdir(fullfile(project_path_full, folders_to_delete{i}), 's');
                                end
                                mkdir(fullfile(project_path_full, folders_to_delete{i}));
                            end

                            % Wipe master folder internals and recreate
                            if exist(fullfile(project_path_full, 'master'), 'dir')
                                rmdir(fullfile(project_path_full, 'master'), 's');
                            end
                            mkdir(fullfile(project_path_full, 'master'));

                            % 2. CREATE project_master.conf
                            projectfolder = strcat('PROJECTFOLDER=',project_path_full);
                            graphsfolder = strcat('GRAPHSFOLDER=', app.InstallRoot, par, 'PHASE_Preprocessing', par, 'snap2stamps', par, 'graphs');

                            if ispc
                                f_project_conf_master = fopen(strcat(project_path_full, '\snap2stamps\bin\project_master.conf'),'w');
                                j = 0;
                                while f_project_conf_master == -1 && j < 20
                                    f_project_conf_master = fopen(strcat(project_path_full, '\snap2stamps\bin\project_master.conf'),'w');
                                    j = j+1;
                                end
                                fprintf(f_project_conf_master,'%s\r\n', '######### CONFIGURATION FILE ######');
                                fprintf(f_project_conf_master,'%s\r\n', '# PROJECT DEFINITION');
                                fprintf(f_project_conf_master,'%s\r\n', projectfolder);
                                fprintf(f_project_conf_master,'%s\r\n', graphsfolder);
                                fprintf(f_project_conf_master,'%s\r\n', '# SPLIT BBOX DEFINITION');
                                fprintf(f_project_conf_master,'LONMIN=%.3f\r\n', lon_min);
                                fprintf(f_project_conf_master,'LATMIN=%.3f\r\n', lat_min);
                                fprintf(f_project_conf_master,'LONMAX=%.3f\r\n', lon_max);
                                fprintf(f_project_conf_master,'LATMAX=%.3f\r\n', lat_max);
                                fprintf(f_project_conf_master,'AUTO_MASTER=%.0f\r\n', auto_master);
                                fprintf(f_project_conf_master,'MASTER_DATE=%s\r\n', master_date);
                                fprintf(f_project_conf_master,'%s\r\n', '# SNAP GPT');
                                fprintf(f_project_conf_master,'GPTBIN_PATH=%s\r\n', gptbin_path);
                                fprintf(f_project_conf_master,'%s\r\n', '# COMPUTING RESOURCES');
                                fprintf(f_project_conf_master,'CPU=%.0f\r\n', cpu);
                                fprintf(f_project_conf_master,'CACHE=%s\r\n', cache);
                                fclose(f_project_conf_master);

                            elseif isunix && ~ismac
                                f_project_conf_master = fopen(strcat(project_path_full, '/snap2stamps/bin/project_master.conf'),'w');
                                j = 0;
                                while f_project_conf_master == -1 && j < 20
                                    f_project_conf_master = fopen(strcat(project_path_full, '/snap2stamps/bin/project_master.conf'),'w');
                                    j = j+1;
                                end
                                fprintf(f_project_conf_master,'%s\n', '######### CONFIGURATION FILE ######');
                                fprintf(f_project_conf_master,'%s\n', '# PROJECT DEFINITION');
                                fprintf(f_project_conf_master,'%s\n', projectfolder);
                                fprintf(f_project_conf_master,'%s\n', graphsfolder);
                                fprintf(f_project_conf_master,'%s\n', '# SPLIT BBOX DEFINITION');
                                fprintf(f_project_conf_master,'LONMIN=%.3f\n', lon_min);
                                fprintf(f_project_conf_master,'LATMIN=%.3f\n', lat_min);
                                fprintf(f_project_conf_master,'LONMAX=%.3f\n', lon_max);
                                fprintf(f_project_conf_master,'LATMAX=%.3f\n', lat_max);
                                fprintf(f_project_conf_master,'AUTO_MASTER=%.0f\n', auto_master);
                                fprintf(f_project_conf_master,'MASTER_DATE=%s\n', master_date);
                                fprintf(f_project_conf_master,'%s\n', '# SNAP GPT');
                                fprintf(f_project_conf_master,'GPTBIN_PATH=%s\n', gptbin_path);
                                fprintf(f_project_conf_master,'%s\n', '# COMPUTING RESOURCES');
                                fprintf(f_project_conf_master,'CPU=%.0f\n', cpu);
                                fprintf(f_project_conf_master,'CACHE=%s\n', cache);
                                fclose(f_project_conf_master);

                            elseif ismac
                                f_project_conf_master = fopen(strcat(project_path_full, '/snap2stamps/bin/project_master.conf'),'w');
                                j = 0;
                                while f_project_conf_master == -1 && j < 20
                                    f_project_conf_master = fopen(strcat(project_path_full, '/snap2stamps/bin/project_master.conf'),'w');
                                    j = j+1;
                                end
                                fprintf(f_project_conf_master,'%s\n', '######### CONFIGURATION FILE ######');
                                fprintf(f_project_conf_master,'%s\n', '# PROJECT DEFINITION');
                                fprintf(f_project_conf_master,'%s\n', projectfolder);
                                fprintf(f_project_conf_master,'%s\n', graphsfolder);
                                fprintf(f_project_conf_master,'%s\n', '# SPLIT BBOX DEFINITION');
                                fprintf(f_project_conf_master,'LONMIN=%.3f\n', lon_min);
                                fprintf(f_project_conf_master,'LATMIN=%.3f\n', lat_min);
                                fprintf(f_project_conf_master,'LONMAX=%.3f\n', lon_max);
                                fprintf(f_project_conf_master,'LATMAX=%.3f\n', lat_max);
                                fprintf(f_project_conf_master,'AUTO_MASTER=%.0f\n', auto_master);
                                fprintf(f_project_conf_master,'MASTER_DATE=%s\n', master_date);
                                fprintf(f_project_conf_master,'%s\n', '# SNAP GPT');
                                fprintf(f_project_conf_master,'GPTBIN_PATH=%s\n', gptbin_path);
                                fprintf(f_project_conf_master,'%s\n', '# COMPUTING RESOURCES');
                                fprintf(f_project_conf_master,'CPU=%.0f\n', cpu);
                                fprintf(f_project_conf_master,'CACHE=%s\n', cache);
                                fclose(f_project_conf_master);
                            end

                            % 3. EXECUTE MASTER SELECTOR AND MASTER PREP BATCH SCRIPT
                            updateOutput(app, 'Running automated CSK Master Auto-Selector and Subsetting...');
                            step_selector = phase_preprocessing_beta.scriptCommand(python, 'CSK_master_selector.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project_master.conf'), app.InstallRoot);
                            step_selector_cmd = step_selector;
                            step1_master = phase_preprocessing_beta.scriptCommand(python, 'CSK_subset_master.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project_master.conf'), app.InstallRoot);
                            step_master_2 = step1_master;

                            path_1_master = fullfile(project_path_full, par, 'snap2stamps', par, 'bin');

                            if ispc
                                f_snap2stamps_master = fopen(strcat(project_path_full, '\snap2stamps\bin\snap2stamps_master.bat'),'w');
                                j = 0;
                                while f_snap2stamps_master == -1 && j < 20
                                    f_snap2stamps_master = fopen(strcat(project_path_full, '\snap2stamps\bin\snap2stamps_master.bat'),'w');
                                    j = j+1;
                                end
                                fprintf(f_snap2stamps_master,'@echo off \r\n');
                                fprintf(f_snap2stamps_master,'cd "%s"\r\n',path_1_master);
                                fprintf(f_snap2stamps_master,'%s\r\n',step_selector_cmd);
                                fprintf(f_snap2stamps_master,'%s\r\n',step_master_2);
                                fclose(f_snap2stamps_master);

                                path_2_master = fullfile(project_path_full, '\snap2stamps\bin\snap2stamps_master.bat');
                                phase_preprocessing_beta.runCommandLive(app, path_2_master, ...
                                    'Master selection and preparation', 3, 18, 1, 1);

                            elseif isunix && ~ismac
                                f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                j = 0;
                                while f_snap2stamps_master == -1 && j < 20
                                    f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                    j = j+1;
                                end
                                fprintf(f_snap2stamps_master,'#!/bin/bash \nset -e\n');
                                fprintf(f_snap2stamps_master,'cd "%s"\n',path_1_master);
                                fprintf(f_snap2stamps_master,'%s\n',step_selector_cmd);
                                fprintf(f_snap2stamps_master,'%s\n',step_master_2);
                                fclose(f_snap2stamps_master);

                                path_2_master = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'));
                                chmod = ['chmod +x' space path_2_master];
                                phase_preprocessing_beta.runCommandLive(app, path_2_master, ...
                                    'Master selection and preparation', 3, 18, 1, 1);

                            elseif ismac
                                f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                j = 0;
                                while f_snap2stamps_master == -1 && j < 20
                                    f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                    j = j+1;
                                end
                                fprintf(f_snap2stamps_master,'#!/bin/bash \nset -e\n');
                                fprintf(f_snap2stamps_master,'cd "%s"\n',path_1_master);
                                fprintf(f_snap2stamps_master,'%s\n',step_selector_cmd);
                                fprintf(f_snap2stamps_master,'%s\n',step_master_2);
                                fclose(f_snap2stamps_master);

                                path_2_master = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'));
                                chmod = ['chmod +x' space path_2_master];
                                phase_preprocessing_beta.runCommandLive(app, path_2_master, ...
                                    'Master selection and preparation', 3, 18, 1, 1);
                            end

                            % FLATTEN DIRECTORY: Move the .h5 file up one level to clean the structure
                            master_subdirs = dir(fullfile(project_path_full, 'master'));
                            for idx = 1:length(master_subdirs)
                                % Check for valid subdirectories, specifically EXCLUDING '.data' folders
                                if master_subdirs(idx).isdir && ~strcmp(master_subdirs(idx).name, '.') && ~strcmp(master_subdirs(idx).name, '..') && ~endsWith(master_subdirs(idx).name, '.data')
                                    subfolder_path = fullfile(master_subdirs(idx).folder, master_subdirs(idx).name);

                                    % Move .h5 up to the root of /master
                                    h5_contents = dir(fullfile(subfolder_path, '*.h5'));
                                    for j = 1:length(h5_contents)
                                        movefile(fullfile(h5_contents(j).folder, h5_contents(j).name), fullfile(project_path_full, 'master'));
                                    end

                                    % Delete the now-empty date folder
                                    rmdir(subfolder_path, 's');
                                end
                            end

                            % 4. VALIDATE OUTPUT
                            master_dim_struct = dir(fullfile(project_path_full, 'master', '*_sub.dim'));
                            if isempty(master_dim_struct)
                                msg = 'ERROR: the master image has not been correctly processed. Check the snap2stamps logs.';
                                error(msg);
                            end
                            master_file_destination = fullfile(master_dim_struct(1).folder, master_dim_struct(1).name);

                            updateOutput(app, '----------------------- STEP 1: Master processing finished -----------------------');

                        else
                            % IF MASTER PROCESSING IS SKIPPED, DYNAMICALLY FIND THE EXISTING MASTER
                            master_dim_struct = dir(fullfile(project_path_full, 'master', '*_sub.dim'));
                            if isempty(master_dim_struct)
                                msg = 'ERROR: No processed master found in /master folder! You must run Master Processing first.';
                                error(msg);
                            end
                            master_file_destination = fullfile(master_dim_struct(1).folder, master_dim_struct(1).name);

                            updateOutput(app, '----------------------- STEP 1: Master processing skipped -----------------------');

                        end

                        if app.StopFlag
                            break;
                        end


                        %% ------------------------ SLAVES PROCESSING -----------------------------

                        % CREATE THE PROJECT.CONF FILE FOR SNAP2STAMPS SLAVES PROCESSING
                        projectfolder = strcat('PROJECTFOLDER=',project_path_full); % path of the project folder
                        graphsfolder = strcat('GRAPHSFOLDER=', app.InstallRoot, par, 'PHASE_Preprocessing', par, 'snap2stamps', par, 'graphs'); % path of the graphs folder of snap2stamps
                        masterfolder = strcat('MASTER=', master_file_destination); % master image path

                        if ispc
                            f_project_conf_slaves = fopen(strcat(project_path_full, '\snap2stamps\bin\project.conf'),'w');
                            j = 0;
                            while f_project_conf_slaves == -1 && j < 20
                                f_project_conf_slaves = fopen(strcat(project_path_full, '\snap2stamps\bin\project.conf'),'w');
                                j = j+1;
                            end
                            fprintf(f_project_conf_slaves,'%s\r\n', '######### CONFIGURATION FILE ######');
                            fprintf(f_project_conf_slaves,'%s\r\n', '###################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# PROJECT DEFINITION');
                            fprintf(f_project_conf_slaves,'%s\r\n', projectfolder);
                            fprintf(f_project_conf_slaves,'%s\r\n', graphsfolder);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# PROCESSING PARAMETERS');
                            fprintf(f_project_conf_slaves,'%s\r\n', masterfolder);
                            fprintf(f_project_conf_slaves,'TC_COHERENCE=%.0f\n', coherence_tc);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# DEM DEFINITION');
                            fprintf(f_project_conf_slaves,'DEMNAME=%s\r\n', dem_name);
                            fprintf(f_project_conf_slaves,'DEMFILE=%s\r\n', dem_file);
                            fprintf(f_project_conf_slaves,'NUMGCP=%.0f\n', num_gcp);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# SPLIT BBOX DEFINITION');
                            fprintf(f_project_conf_slaves,'LONMIN=%.3f\r\n', lon_min);
                            fprintf(f_project_conf_slaves,'LATMIN=%.3f\r\n', lat_min);
                            fprintf(f_project_conf_slaves,'LONMAX=%.3f\r\n', lon_max);
                            fprintf(f_project_conf_slaves,'LATMAX=%.3f\r\n', lat_max);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# SNAP GPT');
                            fprintf(f_project_conf_slaves,'GPTBIN_PATH=%s\r\n', gptbin_path);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# COMPUTING RESOURCES TO EMPLOY');
                            fprintf(f_project_conf_slaves,'CPU=%.0f\r\n', cpu);
                            fprintf(f_project_conf_slaves,'CACHE=%s\r\n', cache);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            f_project_conf_slaves = fclose(f_project_conf_slaves);
                        elseif isunix && ~ismac
                            f_project_conf_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/project.conf'),'w');
                            j = 0;
                            while f_project_conf_slaves == -1 && j < 20
                                f_project_conf_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/project.conf'),'w');
                                j = j+1;
                            end
                            fprintf(f_project_conf_slaves,'%s\r\n', '######### CONFIGURATION FILE ######');
                            fprintf(f_project_conf_slaves,'%s\r\n', '###################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# PROJECT DEFINITION');
                            fprintf(f_project_conf_slaves,'%s\r\n', projectfolder);
                            fprintf(f_project_conf_slaves,'%s\r\n', graphsfolder);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# PROCESSING PARAMETERS');
                            fprintf(f_project_conf_slaves,'%s\r\n', masterfolder);
                            fprintf(f_project_conf_slaves,'TC_COHERENCE=%.0f\n', coherence_tc);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# DEM DEFINITION');
                            fprintf(f_project_conf_slaves,'DEMNAME=%s\r\n', dem_name);
                            fprintf(f_project_conf_slaves,'DEMFILE=%s\r\n', dem_file);
                            fprintf(f_project_conf_slaves,'NUMGCP=%.0f\n', num_gcp);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# SPLIT BBOX DEFINITION');
                            fprintf(f_project_conf_slaves,'LONMIN=%.3f\r\n', lon_min);
                            fprintf(f_project_conf_slaves,'LATMIN=%.3f\r\n', lat_min);
                            fprintf(f_project_conf_slaves,'LONMAX=%.3f\r\n', lon_max);
                            fprintf(f_project_conf_slaves,'LATMAX=%.3f\r\n', lat_max);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# SNAP GPT');
                            fprintf(f_project_conf_slaves,'GPTBIN_PATH=%s\r\n', gptbin_path);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# COMPUTING RESOURCES TO EMPLOY');
                            fprintf(f_project_conf_slaves,'CPU=%.0f\r\n', cpu);
                            fprintf(f_project_conf_slaves,'CACHE=%s\r\n', cache);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            f_project_conf_slaves = fclose(f_project_conf_slaves);
                        elseif ismac
                            f_project_conf_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/project.conf'),'w');
                            j = 0;
                            while f_project_conf_slaves == -1 && j < 20
                                f_project_conf_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/project.conf'),'w');
                                j = j+1;
                            end
                            fprintf(f_project_conf_slaves,'%s\r\n', '######### CONFIGURATION FILE ######');
                            fprintf(f_project_conf_slaves,'%s\r\n', '###################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# PROJECT DEFINITION');
                            fprintf(f_project_conf_slaves,'%s\r\n', projectfolder);
                            fprintf(f_project_conf_slaves,'%s\r\n', graphsfolder);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# PROCESSING PARAMETERS');
                            fprintf(f_project_conf_slaves,'%s\r\n', masterfolder);
                            fprintf(f_project_conf_slaves,'TC_COHERENCE=%.0f\n', coherence_tc);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# DEM DEFINITION');
                            fprintf(f_project_conf_slaves,'DEMNAME=%s\r\n', dem_name);
                            fprintf(f_project_conf_slaves,'DEMFILE=%s\r\n', dem_file);
                            fprintf(f_project_conf_slaves,'NUMGCP=%.0f\n', num_gcp);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# SPLIT BBOX DEFINITION');
                            fprintf(f_project_conf_slaves,'LONMIN=%.3f\r\n', lon_min);
                            fprintf(f_project_conf_slaves,'LATMIN=%.3f\r\n', lat_min);
                            fprintf(f_project_conf_slaves,'LONMAX=%.3f\r\n', lon_max);
                            fprintf(f_project_conf_slaves,'LATMAX=%.3f\r\n', lat_max);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# SNAP GPT');
                            fprintf(f_project_conf_slaves,'GPTBIN_PATH=%s\r\n', gptbin_path);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            fprintf(f_project_conf_slaves,'%s\r\n', '# COMPUTING RESOURCES TO EMPLOY');
                            fprintf(f_project_conf_slaves,'CPU=%.0f\r\n', cpu);
                            fprintf(f_project_conf_slaves,'CACHE=%s\r\n', cache);
                            fprintf(f_project_conf_slaves,'%s\r\n', '##################################');
                            f_project_conf_slaves = fclose(f_project_conf_slaves);
                        end

                        % EXECUTE THE .CONF FILE VIA A BATCH/BASH FILE
                        dp = phase_preprocessing_beta.skipPrefix();
                        step1_slaves = phase_preprocessing_beta.scriptCommand(python, 'CSK_slaves_prep.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % slaves preparation
                        step2_slaves = phase_preprocessing_beta.scriptCommand(python, 'CSK_subset_slaves.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % slaves splitting & apply orbits
                        step3_slaves = phase_preprocessing_beta.scriptCommand(python, 'CSK_coreg_ifg.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % coregistration & interferogram
                        step4_slaves = phase_preprocessing_beta.scriptCommand(python, 'CSK_stamps_export.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % StaMPS export
                        step5_slaves = phase_preprocessing_beta.scriptCommand(python, 'CSK_average_intensity.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % terrain corrected coherence and lia
                        step6_slaves = phase_preprocessing_beta.scriptCommand(python, 'CSK_terrain_correction.py', fullfile(project_path_full, 'snap2stamps', 'bin', 'project.conf'), app.InstallRoot); % terrain corrected coherence and lia

                            % CASES DEFENDING ON FIRST STEP
                            if isnumeric(first_step)
                                first_step_num = first_step;
                            elseif ischar(first_step)
                                first_step_num = str2double(first_step);
                            end
                            if first_step_num == 1
                                step_slaves_1 = step1_slaves;
                                step_slaves_2 = step2_slaves;
                                step_slaves_3 = step3_slaves;
                                step_slaves_4 = step4_slaves;
                                step_slaves_5 = step5_slaves;
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                end
                            elseif first_step_num == 2
                                step_slaves_1 = [dp step1_slaves];
                                step_slaves_2 = step2_slaves;
                                step_slaves_3 = step3_slaves;
                                step_slaves_4 = step4_slaves;
                                step_slaves_5 = step5_slaves;
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                end
                            elseif first_step_num == 3
                                step_slaves_1 = [dp step1_slaves];
                                step_slaves_2 = [dp step2_slaves];
                                step_slaves_3 = step3_slaves;
                                step_slaves_4 = step4_slaves;
                                step_slaves_5 = step5_slaves;
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                end
                            elseif first_step_num == 4
                                step_slaves_1 = [dp step1_slaves];
                                step_slaves_2 = [dp step2_slaves];
                                step_slaves_3 = [dp step3_slaves];
                                step_slaves_4 = step4_slaves;
                                step_slaves_5 = step5_slaves;
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                end
                            elseif first_step_num == 5
                                step_slaves_1 = [dp step1_slaves];
                                step_slaves_2 = [dp step2_slaves];
                                step_slaves_3 = [dp step3_slaves];
                                step_slaves_4 = [dp step4_slaves];
                                step_slaves_5 = step5_slaves;
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                end
                            elseif first_step_num == 6
                                step_slaves_1 = [dp step1_slaves];
                                step_slaves_2 = [dp step2_slaves];
                                step_slaves_3 = [dp step3_slaves];
                                step_slaves_4 = [dp step4_slaves];
                                step_slaves_5 = [dp step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = step6_slaves;
                                else
                                    updateOutput(app, ['You are running just the terrain correction for coherence and LIA bands, but you have set tc_coherence = 1. ' ...
                                        'Please change it to 0 to perform this step.'])
                                end
                            end

                        if update_processed_run
                            step_slaves_5 = [dp step5_slaves];
                            updateOutput(app, 'Update mode: stack-wide average intensity deferred until existing products are restored.');
                        end

                        path_1_slaves = fullfile(project_path_full, par, 'snap2stamps', par, 'bin');
                        slaves_proof = strcat(project_path_full, par, 'slaves', par, 'dummyDownload.txt');

                        if ispc
                            f_snap2stamps_slaves = fopen(strcat(project_path_full, '\snap2stamps\bin\snap2stamps_slaves.bat'),'w');
                            j = 0;
                            while f_snap2stamps_slaves == -1 && j < 20
                                f_snap2stamps_slaves = fopen(strcat(project_path_full, '\snap2stamps\bin\snap2stamps_slaves.bat'),'w');
                                j = j+1;
                            end
                            fprintf(f_snap2stamps_slaves,'@echo off \r\n');
                            fprintf(f_snap2stamps_slaves,'cd "%s"\r\n',path_1_slaves);
                            fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_1);
                            fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_2);
                            fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_3);
                            fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_4);
                            fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_5);
                            if coherence_tc == 0
                                fprintf(f_snap2stamps_slaves,'%s\r\n',step_slaves_6);
                            end
                            fprintf(f_snap2stamps_slaves,'cd "%s"\r\n',path_1_download);
                            fprintf(f_snap2stamps_slaves,'type nul > dummyDownload.txt \r\n');
                            fprintf(f_snap2stamps_slaves,'exit \r\n');
                            f_snap2stamps_slaves = fclose(f_snap2stamps_slaves);
                        elseif isunix && ~ismac
                            f_snap2stamps_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh'),'w');
                            j = 0;
                            while f_snap2stamps_slaves == -1 && j < 20
                                f_snap2stamps_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh'),'w');
                                j = j+1;
                            end
                            fprintf(f_snap2stamps_slaves,'#!/bin/bash \nset -e\n');
                            fprintf(f_snap2stamps_slaves,'cd "%s"\n',path_1_slaves);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_1);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_2);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_3);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_4);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_5);
                            if coherence_tc == 0
                                fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_6);
                            end
                            fprintf(f_snap2stamps_slaves,'sleep 5 \n');
                            fprintf(f_snap2stamps_slaves,'cd "%s"\n',path_1_download);
                            fprintf(f_snap2stamps_slaves,'touch dummyDownload.txt \n');
                            f_snap2stamps_slaves = fclose(f_snap2stamps_slaves);
                        elseif ismac
                            f_snap2stamps_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh'),'w');
                            j = 0;
                            while f_snap2stamps_slaves == -1 && j < 20
                                f_snap2stamps_slaves = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh'),'w');
                                j = j+1;
                            end
                            fprintf(f_snap2stamps_slaves,'#!/bin/bash \nset -e\n');
                            fprintf(f_snap2stamps_slaves,'cd "%s"\n',path_1_slaves);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_1);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_2);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_3);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_4);
                            fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_5);
                            if coherence_tc == 0
                                fprintf(f_snap2stamps_slaves,'%s\n',step_slaves_6);
                            end
                            fprintf(f_snap2stamps_slaves,'sleep 5 \n');
                            fprintf(f_snap2stamps_slaves,'cd "%s"\n',path_1_download);
                            fprintf(f_snap2stamps_slaves,'touch dummyDownload.txt \n');
                            f_snap2stamps_slaves = fclose(f_snap2stamps_slaves);
                        end

                        updateOutput(app, '----------------------- STEP 2: Slaves processing started -----------------------');

                        if ispc
                            path_2_slaves = fullfile(project_path_full, '\snap2stamps\bin\snap2stamps_slaves.bat');
                            phase_preprocessing_beta.runCommandLive(app, path_2_slaves, ...
                                'Slave processing pipeline', 18, 92, first_step_num, 6);
                        elseif isunix && ~ismac
                            path_2_slaves = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh'));
                            chmod_in = ('chmod +x');
                            chmod = [chmod_in space path_2_slaves];
                            xterm = ('sudo xterm -e');

                            phase_preprocessing_beta.runCommandLive(app, path_2_slaves, ...
                                'Slave processing pipeline', 18, 92, first_step_num, 6);
                        elseif ismac
                            path_2_slaves = strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh');
                            chmod_in = 'chmod +x';
                            chmod = [chmod_in space path_2_slaves];
                            % Executed without opening Terminal by runCommandLive.

                            phase_preprocessing_beta.runCommandLive(app, path_2_slaves, ...
                                'Slave processing pipeline', 18, 92, first_step_num, 6);
                        end

                        % wait until the cmd has finished the work
                        while exist(slaves_proof,'file') == 0
                              pause(1);
                        end
                        delete(slaves_proof);

                        % DYNAMIC WAIT FOR GEO-TIFF CREATION (Timeout: 5 minutes)
                        %updateOutput(app, 'Waiting for .tif files to be written on disk...');
                        timeout = 300; % Maximum wait time in seconds
                        elapsed = 0;
                        intensity_folder = fullfile(project_path_full, 'intensity');

                        % Wait until the folder exists and is not empty, but do not exceed the timeout
                        while (~exist(intensity_folder, 'dir') || isempty(dir(fullfile(intensity_folder, '*.tif')))) && elapsed < timeout
                            pause(5); % Check every 5 seconds
                            elapsed = elapsed + 5;
                        end

                        if elapsed >= timeout
                            updateOutput(app, 'WARNING: Timeout exceeded. The .tif files were not found.');
                        else
                            %updateOutput(app, '.tif files found, starting reprojection...');
                            % Leave a 2-second safety margin for the OS to release the file lock
                            pause(2);
                        end

                        % APPLY THE CORRECT PROJECTION TO THE AVERAGE INTENSITY GEOTIFF
                        if first_step_num <= 5
                            intensity_folder = fullfile(project_path_full, 'intensity');
                            if exist(intensity_folder, 'dir')
                                intensity_files = dir(fullfile(intensity_folder, '*.tif'));
                                for i = 1:length(intensity_files)
                                    file_path = fullfile(intensity_folder, intensity_files(i).name);
                                    app.reprojectGeoTiffIfGeographic(file_path, epsg_code);
                                end
                            end
                        end

                        coherence_folder = fullfile(project_path_full, 'coherence');
                        if exist(coherence_folder, 'dir')
                            coherence_files = dir(fullfile(coherence_folder, '*.tif'));
                            for i = 1:length(coherence_files)
                                file_path = fullfile(coherence_folder, coherence_files(i).name);
                                app.reprojectGeoTiffIfGeographic(file_path, epsg_code);
                            end
                        end

                        lia_folder = fullfile(project_path_full, 'lia');
                        if exist(lia_folder, 'dir')
                            lia_files = dir(fullfile(lia_folder, '*.tif'));
                            for i = 1:length(lia_files)
                                file_path = fullfile(lia_folder, lia_files(i).name);
                                app.reprojectGeoTiffIfGeographic(file_path, epsg_code, 2);
                            end
                        end

                        updateOutput(app, '----------------------- STEP 2: Slaves processing finished -----------------------');

                        if app.StopFlag
                            break;
                        end

                        % DYNAMIC STAMPS FOLDER AUTOMATION (ORBIT & DATE EXTRACTION)
                        coreg_dir = fullfile(project_path_full, 'coreg');
                        dim_files = dir(fullfile(coreg_dir, '*.dim'));

                        if isempty(dim_files)
                            updateOutput(app, 'ERROR: No coregistered files found to extract dates/orbit!');
                            return;
                        end

                        % 1. Extract all dates from the filenames to find absolute Min and Max
                        all_dates = NaT(0); % Initialize empty datetime array
                        for idx = 1:length(dim_files)
                            % Parse the YYYYMMDD_YYYYMMDD.dim filename structure
                            tokens = regexp(dim_files(idx).name, '(\d{8})_(\d{8})\.dim', 'tokens');
                            if ~isempty(tokens)
                                all_dates(end+1, 1) = datetime(tokens{1}{1}, 'InputFormat', 'yyyyMMdd');
                                all_dates(end+1, 1) = datetime(tokens{1}{2}, 'InputFormat', 'yyyyMMdd');
                            end
                        end

                        % Sort and format to English MmmYY (e.g., Jul19, Aug23)
                        min_date = min(all_dates);
                        max_date = max(all_dates);
                        date_str_start = char(datetime(min_date, 'Format', 'MMMyy', 'Locale', 'en_US'));
                        date_str_end   = char(datetime(max_date, 'Format', 'MMMyy', 'Locale', 'en_US'));

                        % 2. Read the Ground Truth Orbit from the Master .dim XML
                        dim_path = fullfile(dim_files(1).folder, dim_files(1).name);
                        dim_text = fileread(dim_path);

                        orbit_type = 'UNKNOWN';
                        if contains(upper(dim_text), 'ASCENDING')
                            orbit_type = 'ASC';
                        elseif contains(upper(dim_text), 'DESCENDING')
                            orbit_type = 'DES';
                        end

                        % 3. Automatically Construct the Target Folder Name!
                        stamps_folder = sprintf('%s_%s_%s', orbit_type, date_str_start, date_str_end);
                        updateOutput(app, ['Dynamically generated StaMPS folder name: ' stamps_folder]);

                        project_parent_path_full = phase_preprocessing_beta.stampsFolder(app.ProjectRoot);
                        stamps_folder_full = fullfile(project_parent_path_full, stamps_folder);
                        if ~isfolder(stamps_folder_full)
                            mkdir(stamps_folder_full);
                            updateOutput(app, ['Created StaMPS dataset folder: ' stamps_folder_full]);
                        end
                        % The beta launches the editable StaMPS module from its
                        % canonical installation. No PHASE_StaMPS.mlapp is copied
                        % or required at runtime.
                        stamps_diagnostic = fullfile(app.InstallRoot, 'PHASE_Preprocessing', 'diagnose_PHASE_StaMPS.m');
                        if isfile(stamps_diagnostic)
                            copyfile(stamps_diagnostic, stamps_folder_full, 'f');
                        end

                        % EVENTUAL REMOVAL OF THE SLAVES IMAGES TO SAVE DISK SPACE
                        if slaves_removal == 0 % check if the removal of the slaves images has to be done or not

                            updateOutput(app, '----------------------- STEP 3: Slaves removal started -----------------------');

                            rmdir(strcat(project_path_full, par, 'slaves'), 's'); % removal of slaves folder to save disk space
                            mkdir(project_path_full, 'slaves'); % create slaves folder

                            updateOutput(app, '----------------------- STEP 3: Slaves removal finished -----------------------');

                            else

                            updateOutput(app, '----------------------- STEP 3: Slaves removal skipped -----------------------');

                        end

                        % Seed a new dataset with the installed configuration, but never
                        % overwrite an existing input_StaMPS.mat: it may contain parameters
                        % the user tuned specifically for this processing run.
                        src_input_mat = fullfile(project_path_full, 'input_StaMPS.mat');
                        dst_input_mat = fullfile(stamps_folder_full, 'input_StaMPS.mat');
                        if isfile(dst_input_mat)
                            updateOutput(app, 'Preserved the existing dataset input_StaMPS.mat.');
                        elseif isfile(src_input_mat)
                            copyfile(src_input_mat, dst_input_mat);
                            updateOutput(app, 'Copied input_StaMPS.mat to the new StaMPS processing folder.');
                        else
                            updateOutput(app, ['NOTICE: input_StaMPS.mat is not available yet. ', ...
                                'PHASE_StaMPS_beta will create it when the initial settings are saved.']);
                        end

                        stamps_app_full = fullfile(project_parent_path_full, stamps_folder);
                        updateOutput(app, ['Preprocessing completed. StaMPS dataset folder: ' stamps_app_full]);
                        updateOutput(app, 'Open the StaMPS PSI section in the PHASE hub when ready.');

                        %% ----------------------------------------------------

                        % Display a message when the script is done
                        updateOutput(app, 'Script CSK_Preprocessing.m completed.');

                        % Change lamp color to green (success)

                    catch ME

                        % If an error occurs
                        if strcmp(ME.identifier,'PHASE:ProcessingStopped')
                            updateOutput(app, 'Preprocessing was force-stopped by the user.');
                        else
                            updateOutput(app, 'An error occurred during script execution. Please check each step log file!');
                        end

                        % Change lamp color to red (error)

                        % Show the actual error
                        rethrow(ME)

                    end

                    condition = 1;

                end

            end

        end

    end
end

function name = shortFileName(pathValue)
[~,stem,extension] = fileparts(pathValue);
name = [stem extension];
end
