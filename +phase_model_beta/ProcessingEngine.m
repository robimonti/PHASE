classdef ProcessingEngine < dynamicprops
    %PROCESSINGENGINE Native .m-only backend for PHASE displacement modelling.
    % The visible App owns the UI; this object owns configuration, progress,
    % cancellation and the unchanged scientific modelling sequence.

    properties
        ConfigRoot = ''
        OutputRoot = ''
        outputDir = ''
        StopRequested = false
        ExternalLogCallback = []
        ExternalProgressCallback = []
    end

    methods
        function app = ProcessingEngine(config,configRoot,outputRoot)
            if nargin < 1 || isempty(config)
                config = phase_model_beta.defaultConfig();
            end
            if nargin >= 2, app.ConfigRoot = char(string(configRoot)); end
            if nargin >= 3, app.OutputRoot = char(string(outputRoot)); end
            app.applyConfig(config);
        end

        function applyConfig(app,config)
            if ~isstruct(config) || ~isscalar(config)
                error('PHASE_Model_beta:invalidEngineConfig', ...
                    'The model configuration must be a scalar structure.');
            end
            names = fieldnames(config);
            for k = 1:numel(names)
                name = names{k};
                if ~isprop(app,name), addprop(app,name); end
                app.(name) = config.(name);
            end
            app.pythonPath = phase_model_beta.resolvePythonPath(app.pythonPath);
        end

        function notifyBetaProgress(app,percentage,phase)
            progress = struct('percentage',double(percentage), ...
                'phase',char(string(phase)));
            if ~isempty(app.ExternalProgressCallback)
                app.ExternalProgressCallback(progress);
            end
            if ~isempty(app.ExternalLogCallback)
                app.ExternalLogCallback(progress.phase);
            end
            drawnow limitrate
            phase_model_beta.throwIfStopped(app);
        end

        function run(app)


            % set figures always to light mode
            set(groot, 'defaultFigureColor', 'w')     % figure background
            set(groot, 'defaultAxesColor',   'w')     % axes background
            set(groot, 'defaultAxesXColor',  'k')     % axis colors
            set(groot, 'defaultAxesYColor',  'k')
            set(groot, 'defaultAxesZColor',  'k')

            % initialize waitbar
            app.notifyBetaProgress(0, 'Preparing environment...');
            drawnow;

            try

                % Every relative path is rooted explicitly for the
                % complete run. Results are stored beside the visible PHASE
                % shortcuts rather than inside the editable engine clone.
                runtimeRoot = phase_model_beta.projectRoot();
                previousRunFolder = pwd;
                runFolderCleanup = onCleanup(@() cd(previousRunFolder)); %#ok<NASGU>
                cd(runtimeRoot);
                outputRoot = fileparts(runtimeRoot);
                if ~isempty(app.OutputRoot)
                    outputRoot = app.OutputRoot;
                end

                % --- 0. Prepare the environment ---

                % - 0.1) Add the path to the MatlabFunctions folder
                    addpath(fullfile(phase_model_beta.projectRoot(), 'MatlabFunctions'));

                % - 0.2) Detect the current environment
                detectedOS = detectOS();

                % - 0.3) Setup the python environment
                setupPythonEnvironment_app(app.pythonPath);

                % - 0.4) Get the current date
                currentDT = datetime('now');

                % - 0.5) Get the logo filepath
                logo_filename = fullfile(phase_model_beta.projectRoot(), 'PHASE_logo.png');

                % Update progress (10%)
                app.notifyBetaProgress(10, 'Loading inputs...');
                drawnow;

                % --- 1. Required inputs ---
                % - 1.1) Path to the PS displacement time series file
                % It must be the one obtained with PHASE - module 1b
                filepathIN = app.filepathIN;


                % - 1.2) Date of the first displacement
                %        In datetime format
                t0IN = app.t0IN;


                % - 1.3) Flag to remove all existing output folders
                flag_ouputDir = app.flag_ouputDir;


                % - 1.4) AOI for PS selection
                % Needed to distinguish among all PS the ones inside the AOI
                % a) flag for selection between shapefile and manual bounding box
                     % false = shapefile, true = bounding box
                    flag_AOIbb = app.flag_AOIbb;

                    % Initialize ALL AOI variables to empty/NaN first
                    filepathAOI = '';
                    lonMinAOI = NaN; lonMaxAOI = NaN;
                    latMinAOI = NaN; latMaxAOI = NaN;

                % b) retrieve input variables
                switch flag_AOIbb
                    case true
                    lonMinAOI = app.lonMinAOI;
                    lonMaxAOI = app.lonMaxAOI;
                    latMinAOI = app.latMinAOI;
                    latMaxAOI = app.latMaxAOI;
                    case false
                    filepathAOI = app.filepathAOI;
                end


                % - 1.5) Project dimension
                projDim = app.projDim;


                % - 1.6) Processing type
                procType = app.procType;


                % - 1.7) a-priori variance (auto / manual / coherence)
                varNoise_method = app.varNoise_method;

                % Initialize defaults
                varNoise_manual = NaN;
                coherence_dir = '';
                constellation = '';
                num_looks = [];

                    % - auto
                    switch varNoise_method
                    case 'manual'
                    % - manual
                    varNoise_manual = app.varNoise_manual; % [mm^2]
                    case 'coherence'
                    % - coherence
                    coherence_dir = app.coherence_dir;
                    constellation = app.constellation;
                    num_looks = app.num_looks;
                    end


                % - 1.8) number of splines (auto / manual)
                num_spl_method = app.num_spl_method;

                % Initialize defaults
                spline_method = '';
                num_spl_manual = NaN;

                    switch num_spl_method
                    case 'auto'
                    % - auto: variance / MDL / F_test / chi2_test
                    spline_method = app.spline_method;
                    case 'manual'
                    % - manual (minimum is 2)
                    num_spl_manual = app.num_spl_manual;
                    end


                % - 1.9) splines regularization parameter lambda (auto / manual)
                lambda_method = app.lambda_method;
                lambda_manual = NaN;
                    % - auto: implemented in the splines section
                    switch lambda_method
                    case 'manual'
                    % - manual
                    lambda_manual = app.lambda_manual;
                    end


                % - 1.10) choose collocation approach (filtering / prediction)
                coll_proc = app.coll_proc;
                coll_step_est = NaN;
                    switch coll_proc
                        case 'prediction'
                    coll_step_est = app.coll_step_est;
                    end


                % - 1.11) Spatial resolution for 1D interpolation [meters]
                cline_resolution = app.cline_resolution;


                % - 1.12) Spatial resolution for 2D interpolation [meters]
                grid_resolution = app.grid_resolution;


                % - 1.13) Minimum number of months for splines filtering
                minMonths = app.minMonths;


                % - 1.14a) Step to bin temporal covariance in STC1D (NaN means automatic)
                dtCov_STC1D = app.dtCov_STC1D;

                % - 1.14b) Step to bin spatial covariance in STC1D (NaN means automatic)
                dsCov_STC1D = app.dsCov_STC1D;


                % - 1.15) DET1D - a-priori noise
                varNoise_DET1D = app.varNoise_DET1D;       % (auto / manual)
                if strcmp(varNoise_DET1D, 'manual')
                    varNoise_manual_DET1D = app.varNoise_manual_DET1D;
                else
                    varNoise_manual_DET1D = NaN;
                end

                % - 1.16) DET1D - splines
                num_spl_method_DET1D = app.num_spl_method_DET1D;    % (auto / manual)
                spline_method_DET1D = app.spline_method_DET1D;      % (variance / MDL / F_test / chi2_test)
                if strcmp(num_spl_method_DET1D, 'manual')
                    % Read directly from UI
                    num_spl_row_manual_DET1D = app.num_spl_row_manual_DET1D;
                    num_spl_col_manual_DET1D = app.num_spl_col_manual_DET1D;
                else
                    num_spl_row_manual_DET1D = NaN;
                    num_spl_col_manual_DET1D = NaN;
                end

                % - 1.17) DET1D - lambda
                lambda_method_DET1D = app.lambda_method_DET1D;   % (auto / manual)
                if strcmp(lambda_method_DET1D, 'manual')
                     lambda_manual_DET1D = app.lambda_manual_DET1D;
                else
                     lambda_manual_DET1D = NaN;
                end

                % - 1.18) ST - type of spatial means
                useInclinedMeansST = app.useInclinedMeansST;  % true for tilted; false for horizontal

                % - 1.19) Step to bin temporal covariance in STC2D (NaN means automatic)
                dtCov_STC2D = app.dtCov_STC2D;

                % - 1.20) Step to bin spatial covariance in STC2D (NaN means automatic)
                dsCov_STC2D = app.dsCov_STC2D;

                % - 1.21) DET2D - a-priori noise
                varNoise_DET2D = app.varNoise_DET2D;       % (auto / manual)
                if strcmp(varNoise_DET2D, 'manual')
                    varNoise_manual_DET2D = app.varNoise_manual_DET2D;
                else
                    varNoise_manual_DET2D = NaN;
                end

                % - 1.22) DET2D - splines
                num_spl_method_DET2D = app.num_spl_method_DET2D;           % (auto / manual)
                spline_method_DET2D = app.spline_method_DET2D;             % (variance / MDL / F_test / chi2_test)
                if strcmp(num_spl_method_DET2D, 'manual')
                    % Read directly from UI
                    num_spl_row_manual_DET2D = app.num_spl_row_manual_DET2D;
                    num_spl_col_manual_DET2D = app.num_spl_col_manual_DET2D;
                    num_spl_t_manual_DET2D   = app.num_spl_t_manual_DET2D;
                else
                    num_spl_row_manual_DET2D = NaN;
                    num_spl_col_manual_DET2D = NaN;
                    num_spl_t_manual_DET2D   = NaN;
                end

                % - 1.23) DET2D - lambda
                lambda_method_DET2D = app.lambda_method_DET2D;      % (auto / manual)
                if strcmp(lambda_method_DET2D, 'manual')
                    lambda_manual_DET2D = app.lambda_manual_DET2D;
                else
                    lambda_manual_DET2D = NaN;
                end

                % - 1.25) Time estimation step
                step_t_ST = app.step_t_ST;

                % - 1.26) covariance model for time modelling (all STC)
                tCovModel_STC1D = app.tCovModel_STC1D;      % alternative is gaussian_cos / exponential
                tCovModel_STC2D = app.tCovModel_STC2D;      % alternative is gaussian_cos / exponential

                % - 1.27) covariance model for space modelling (all STC)
                sCovModel_STC1D = app.sCovModel_STC1D;      % alternative is gaussian_cos / exponential
                sCovModel_STC2D = app.sCovModel_STC2D;      % alternative is gaussian_cos / exponential

                % - 1.26) Handling for spatial correlation
                detrendMethodST = app.detrendMethod;

                % - 1.28) maximum polynomial degree (all) - options: 1 / 2 / 3
                polyDegreeST = str2double(app.polyDegreeST);

                % - 1.29) marker size for geoplots (all)
                markerSize = app.markerSize;

                % - 1.30) Interpolation of time series
                flag_PSinterp = app.flag_PSinterp;
                flag_tsExtr = app.flag_tsExtr;
                filepath_EXTR = app.filepath_EXTR;

                % Update progress (20%)
                app.notifyBetaProgress(20, 'Importing files...');
                drawnow;



                %% --- 2. File import ---
                % Import the PS displacement time series

                % - 2.1a) .xlsx format
                if contains(filepathIN, '.xlsx', 'IgnoreCase', true)
                    fileIN = readmatrix(filepathIN);
                    dataIN = fileIN(3:end, :);
                    displIN = dataIN(:, 5:end);
                    t_relIN = fileIN(2, 5:end);
                    t_dateIN = t_relIN + t0IN;
                    PSidIN = dataIN(:,1);
                    lonlatIN = dataIN(:,2:3);

                % - 2.1b) .csv format
                elseif contains(filepathIN, '.csv', 'IgnoreCase', true)
                    fileIN = readmatrix(filepathIN, 'DecimalSeparator', ',');
                    dataIN = fileIN(3:end, :);
                    displIN = dataIN(:, 5:end);
                    t_relIN = fileIN(2, 5:end);
                    t_dateIN = t_relIN + t0IN;
                    PSidIN = dataIN(:,1);
                    lonlatIN = dataIN(:,2:3);

                else
                % - 2.1c) Invalid file type
                    error('The file must be either a .xlsx or .csv file.');
                end


                % - 2.2) Import the shapefile of the AOI
                % Loaded later by phase_model_beta.readAoiShapefile so the
                % map and numerical selection share one geometry.


                % - 2.3) Import coordinates of unwrapping reference point
                [parentFolder_tmp, filenameIN, ~] = fileparts(filepathIN);
                parentFolder = fileparts(parentFolder_tmp);
                clear parentFolder_tmp
                filepathREF = fullfile(parentFolder, 'input_StaMPS.mat');

                % check if the .mat file exists
                if isfile(filepathREF)
                    % load variables if the file exists
                    load(filepathREF, 'ref_centre_lonlat', 'ref_radius');

                    % validate reference coordinates
                    if ref_centre_lonlat(1) == 0 && ref_centre_lonlat(2) == 0
                        ref_centre_lonlat = [NaN, NaN];
                        ref_radius = NaN;
                    end
                else
                    % assign NaN if the file does not exist
                    ref_centre_lonlat = [NaN, NaN];
                    ref_radius = NaN;
                end

                % make a circle around the reference point if valid
                if ~any(isnan(ref_centre_lonlat)) && ~isnan(ref_radius)
                    lonlat_circle = plot_circle_geo(ref_centre_lonlat(1), ref_centre_lonlat(2), ref_radius);
                else
                    lonlat_circle = [];
                end

                % Update progress (30%)
                app.notifyBetaProgress(30, 'Creating folders...');
                drawnow;



                %% --- 3. Check/Create output folder ---
                % Folder structure to store processing files

                % - 3.1) Removal of all output folders
                % get the status quo
                baseFolderName = 'output';
                allFolders = dir(outputRoot);
                allFolderNames = {allFolders([allFolders.isdir]).name};
                outputFolders = allFolderNames(startsWith(allFolderNames, baseFolderName));

                % decide for removal
                switch flag_ouputDir
                    case true
                    if ~isempty(outputFolders)
                        fprintf('Removing existing output folders...\n');
                        for i = 1:numel(outputFolders)
                            folderToRemove = outputFolders{i};
                            if ~strcmp(folderToRemove, '.') && ~strcmp(folderToRemove, '..')
                                rmdir(fullfile(outputRoot,folderToRemove), 's');
                                fprintf('Removed folder: %s\n', folderToRemove);
                            end
                        end
                    else
                        fprintf('No existing output folders to remove.\n');
                    end
                    case false
                    fprintf('Flag set to false. Skipping folder removal.\n');
                end


                % - 3.2) Output folder
                % get the status quo
                allFolders = dir(outputRoot);
                allFolderNames = {allFolders([allFolders.isdir]).name};
                outputFolders = allFolderNames(startsWith(allFolderNames, baseFolderName));

                % determine the next folder name
                if isempty(outputFolders)
                    % no "output" folder exists
                    outputDir = sprintf('%s_%03d', baseFolderName, 1);

                else
                    % extract numbers from folder names and find the highest number
                    numbers = cellfun(@(x) sscanf(x, [baseFolderName '_%d']), outputFolders, 'UniformOutput', false);
                    numbers = cell2mat(numbers);
                    if isempty(numbers)
                        % handle edge case if no numbered folders exist
                        outputDir = sprintf('%s_%03d', baseFolderName, 1);
                    else
                        % increment the highest number
                        outputDir = sprintf('%s_%03d', baseFolderName, max(numbers) + 1);
                    end
                end

                % Use the explicit project result root when one is selected.
                outputDir = fullfile(outputRoot,outputDir);
                [created,createMessage] = mkdir(outputDir);
                if ~created
                    error('PHASE_Model_beta:outputCreateFailed', ...
                        'Could not create output folder %s: %s',outputDir,createMessage);
                end
                app.outputDir = char(java.io.File(outputDir).getCanonicalPath());
                fprintf('Output folder created: %s\n',app.outputDir);


                % folder for figures
                figsDir = fullfile(outputDir, 'figures');

                % folder for files
                filesDir = fullfile(outputDir, 'files');


                % - 3.3) Create processing folders for geoSplinter
                gS_input_path = fullfile(outputDir, 'geoSplinter', 'data_input');
                gS_output_path = fullfile(outputDir, 'geoSplinter', 'data_output');
                gS_job_path = fullfile(outputDir, 'geoSplinter', 'job');
                gS_synth_path = fullfile(outputDir, 'geoSplinter', 'data_synthesis');

                processingFolders = {figsDir,filesDir, ...
                    fullfile(filesDir,'shp'),fullfile(filesDir,'mat'), ...
                    gS_input_path,gS_output_path,gS_job_path,gS_synth_path};
                for folderIndex = 1:numel(processingFolders)
                    folder = processingFolders{folderIndex};
                    [folderCreated,folderMessage] = mkdir(folder);
                    if ~folderCreated && ~isfolder(folder)
                        error('PHASE_Model_beta:processingFolderCreateFailed', ...
                            'Could not create modelling folder %s: %s',folder,folderMessage);
                    end
                end


                % Update progress (40%)
                app.notifyBetaProgress(40, 'Preparing data...');
                drawnow;



                %% --- 4. Data preparation ---
                % Visualization of the imported data and prepration of the required
                % variables

                % - 4.1) Geographic AOI, from the map or an explicit .prj
                switch flag_AOIbb
                    case true
                        lonlatAOI = app.aoi_polygon_lonlat;
                        if size(lonlatAOI,2) ~= 2 || size(lonlatAOI,1) < 3
                            lonlatAOI = [
                                lonMinAOI, latMinAOI;
                                lonMaxAOI, latMinAOI;
                                lonMaxAOI, latMaxAOI;
                                lonMinAOI, latMaxAOI;
                                lonMinAOI, latMinAOI
                            ];
                        elseif ~isequal(lonlatAOI(1,:),lonlatAOI(end,:))
                            lonlatAOI(end+1,:) = lonlatAOI(1,:);
                        end
                    case false
                        [lonlatAOI,~,aoiInfo] = phase_model_beta.readAoiShapefile( ...
                            filepathAOI,filepathIN);
                        fprintf('AOI coordinates: %s; polygon parts: %d\n', ...
                            aoiInfo.coordinateType,aoiInfo.partCount);
                end

                % - 4.2) Select the dominant UTM zone once for the AOI.
                [projectCRS,utmZone,epsgCode] = ...
                    phase_model_beta.selectProjectUtm(lonlatAOI);
                [xIN,yIN] = phase_model_beta.projectWgs84(projectCRS,lonlatIN);
                xyIN = [xIN,yIN];
                finiteAOI = all(isfinite(lonlatAOI),2);
                xyAOI = NaN(size(lonlatAOI));
                [xAOI,yAOI] = phase_model_beta.projectWgs84( ...
                    projectCRS,lonlatAOI(finiteAOI,:));
                xyAOI(finiteAOI,:) = [xAOI,yAOI];
                fprintf('Project coordinate reference: EPSG:%d (%s)\n', ...
                    epsgCode,utmZone);

                % check which PS are inside the AOI
                xyIN_AOI_flag = inpolygon(xyIN(:,1), xyIN(:,2), xyAOI(:,1), xyAOI(:,2));
                xyIN_AOI = xyIN(xyIN_AOI_flag, :);
                lonlatIN_AOI = lonlatIN(xyIN_AOI_flag, :);
                dataIN_AOI = dataIN(xyIN_AOI_flag, :);
                displIN_AOI = displIN(xyIN_AOI_flag, :);
                PSidIN_AOI = PSidIN(xyIN_AOI_flag, :);
                if isempty(PSidIN_AOI)
                    error('PHASE_Model_beta:noPsInsideAoi', ...
                        ['The selected AOI contains no persistent scatterers from the ', ...
                         'input dataset. Check the AOI shown on the map or select the ', ...
                         'full PS extent before starting.']);
                end


                % - 4.3) Figure of processing scene & AOI
                f = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
                geobasemap satellite
                hold on
                legendHandles = gobjects(0); legendLabels = {};
                aoiPlot = geoplot(geopolyshape(lonlatAOI(:,2), lonlatAOI(:,1)), ...
                    'FaceColor', '#FFFF9F', 'EdgeColor', 'black', 'LineWidth', 1);
                legendHandles(end+1) = aoiPlot(1);
                legendLabels{end+1} = 'AOI';
                if contains(filepathIN, 'ASC')
                    % detected ASC orbit data
                    legendHandles(end+1) = geoscatter(lonlatIN(:,2), lonlatIN(:,1), 30, 'r', 'filled', 'MarkerEdgeColor', 'k');
                    legendHandles(end+1) = geoscatter(lonlatIN_AOI(:,2), lonlatIN_AOI(:,1), 30, 'm', 'filled', 'MarkerEdgeColor', 'k');
                    title('Imported PS - ASC orbit', 'FontSize', 20)
                elseif contains(filepathIN, 'DSC')
                    % detected DSC orbit data
                    legendHandles(end+1) = geoscatter(lonlatIN(:,2), lonlatIN(:,1), 30, 'Color', [30 144 255]/255, 'MarkerEdgeColor', 'k');
                    legendHandles(end+1) = geoscatter(lonlatIN_AOI(:,2), lonlatIN_AOI(:,1), 30, 'c', 'filled', 'MarkerEdgeColor', 'k');
                    title('Imported PS - DSC orbit', 'FontSize', 20)
                else
                    % no detected orbit
                    legendHandles(end+1) = geoscatter(lonlatIN(:,2), lonlatIN(:,1), 30, 'm');
                    legendHandles(end+1) = geoscatter(lonlatIN_AOI(:,2), lonlatIN_AOI(:,1), 30, 'm', 'filled');
                    title('Imported PS', 'FontSize', 20)
                end
                legendLabels(end+1:end+2) = {'PS outside AOI','PS inside AOI'};
                if ~any(isnan(ref_centre_lonlat)) && ~isnan(ref_radius)
                    % add reference area used for unwrapping
                    legendHandles(end+1) = geoscatter(ref_centre_lonlat(2), ref_centre_lonlat(1), 10, 'g', 'filled');
                    geoplot(lonlat_circle(:,2), lonlat_circle(:,1), 'g', 'LineWidth', 1.2);
                    legendLabels{end+1} = 'unwrapping ref.';
                end
                legend(legendHandles,legendLabels,'FontSize',13)
                fig1_filename = strcat(figsDir, filesep, 'AOI_PS.png');
                phase_model_beta.exportFigure(f,fig1_filename);
                close(f)


                % - 4.4) Create the interpolation grid / centerline based on projDim
                if strcmp(procType, 'temporal')
                    % Pure temporal modelling works at the observed PS only.
                    % Do not allocate an unused centerline or potentially huge 2D grid.
                    centerline_data = []; xy_grid = []; lonlat_grid_AOI = [];
                    fprintf('Pure temporal mode: spatial grid generation skipped.\n');
                else
                if size(xyAOI, 1) < 3
                    error('xyAOI must have at least 3 points to define a polygon.');
                end
                switch projDim
                    case '1D'
                        if cline_resolution <= 0
                            error('cline_resolution must be positive.');
                        end
                    case '2D'
                        if grid_resolution <= 0
                            error('grid_resolution must be positive.');
                        end
                end

                switch projDim
                    case '1D'
                        % Compute centerline using the function
                        [centerline_data, xy_grid, lonlat_grid_AOI] = compute_centerline(xyAOI, cline_resolution, xyIN_AOI);

                    case '2D'
                        % xy grid
                        [x_grid, y_grid] = meshgrid(min(xyAOI(:,1)):grid_resolution:max(xyAOI(:,1)), ...
                                                    min(xyAOI(:,2)):grid_resolution:max(xyAOI(:,2)));

                        % filter points inside AOI
                        inAOI_xy = inpolygon(x_grid(:), y_grid(:), xyAOI(:,1), xyAOI(:,2));
                        x_grid_AOI = x_grid(inAOI_xy);
                        y_grid_AOI = y_grid(inAOI_xy);
                        xy_grid = [x_grid_AOI, y_grid_AOI];

                        % convert to geographic coordinates for plotting
                        [lat_grid_AOI, lon_grid_AOI] = projinv(projectCRS,x_grid_AOI,y_grid_AOI);
                        lonlat_grid_AOI = [lon_grid_AOI, lat_grid_AOI];
                        centerline_data = [];
                end


                end

                % - 4.5) Determine municipality and define export filenames
                % use the reference point if valid; otherwise, use the center of the PS data inside the AOI
                if numel(ref_centre_lonlat)>=2 && all(isfinite(ref_centre_lonlat(1:2)))
                    query_lon = ref_centre_lonlat(1);
                    query_lat = ref_centre_lonlat(2);
                else
                    query_lon = mean(lonlatIN_AOI(:,1), 'omitnan');
                    query_lat = mean(lonlatIN_AOI(:,2), 'omitnan');
                end

                % fetch the location data
                if isfinite(query_lon) && isfinite(query_lat) && ...
                        abs(query_lon)<=180 && abs(query_lat)<=90
                    [municipality, country] = get_place_from_coordinates(query_lon, query_lat);
                else
                    municipality = 'Unknown'; country = 'Unknown';
                    fprintf('Reverse geocoding skipped because the AOI centre is invalid.\n');
                end
                municipality_exp = create_safe_filename(municipality);

                % define the export filename
                if strcmp(municipality_exp, 'untitled') || isempty(municipality_exp)
                    filenameOUT = strcat(filenameIN, '_out');
                else
                    filenameOUT = strcat(municipality_exp, '_out');
                end

                % define the Excel export filename
                filenameOUT_e = fullfile(outputDir, strcat(filenameOUT, '.xlsx'));

                % define the .shp export filename
                filenameOUT_s = fullfile(filesDir, 'shp', strcat(filenameOUT, '.shp'));


                % - 4.6) Figure of AOI context
                % load country shapefile
                countriesShpPath = fullfile('Extra', 'NaturalEarth', 'ne_50m_admin_0_countries.shp');
                countries = shaperead(countriesShpPath, 'UseGeoCoords', true);

                aoi_mean_lon = mean(lonlatAOI(:,1),'omitnan');
                aoi_mean_lat = mean(lonlatAOI(:,2),'omitnan');

                % Pre-allocate to prevent undefined variable crashes
                countryShp = [];

                % 1. Primary Check: Find the exact country containing the AOI
                for i = 1:length(countries)
                    if inpolygon(aoi_mean_lon, aoi_mean_lat, countries(i).Lon, countries(i).Lat)
                        countryShp = countries(i);
                        break;
                    end
                end

                % 2. Fallback: If AOI center is in the ocean/bay, find the closest country border
                if isempty(countryShp)
                    min_dist = inf;
                    closest_idx = 1;
                    for i = 1:length(countries)
                        % Calculate squared distance to the closest border point of each country
                        dist = min((countries(i).Lon - aoi_mean_lon).^2 + (countries(i).Lat - aoi_mean_lat).^2);
                        if dist < min_dist
                            min_dist = dist;
                            closest_idx = i;
                        end
                    end
                    countryShp = countries(closest_idx);
                    fprintf('AOI center fell outside land borders. Snapped to the closest country boundary.\n');
                end

                % plot
                f = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
                geobasemap topographic
                hold on
                geoplot(countryShp.Lat, countryShp.Lon, 'm-', 'LineWidth', 1.5);
                geoplot(aoi_mean_lat, aoi_mean_lon, 'r.', 'MarkerSize', 15);
                text(aoi_mean_lat, aoi_mean_lon, ' AOI', 'FontSize', 20, 'Color', 'red', 'FontWeight', 'bold');
                title('AOI location in its country', 'FontSize', 20)
                fig2_filename = strcat(figsDir, filesep, 'AOI_context.png');
                phase_model_beta.exportFigure(f,fig2_filename);
                close(f)


                % Update progress (50%)
                app.notifyBetaProgress(50, 'Modeling displacement...');
                drawnow;



                %% --- 5. Displacement time series modelling ---
                % Each displacement time series is independently interpolated with a robust
                % deterministic and stochastic methodology. This is chosen when the goal is
                % to have the best interpolation in time, disregarding the spatial position
                % of the PS

                % validate inputs
                if isempty(dataIN_AOI) || isempty(displIN_AOI) || isempty(PSidIN_AOI) || ...
                   isempty(t_dateIN) || isempty(t_relIN) || isempty(xyIN_AOI)
                    error('One or more required inputs for ModellingInTime are empty.');
                end

                % - 5.1) Create OptionalArgs cell array for ModellingInTime function
                    OptionalArgs = {'stop_check', ...
                        @() phase_model_beta.throwIfStopped(app)};
                    switch varNoise_method
                        case 'manual'
                        OptionalArgs = [OptionalArgs, {'varNoise_manual', varNoise_manual}];
                        case 'coherence'
                        OptionalArgs = [OptionalArgs, {'utmZone', utmZone, 'coherence_dir', coherence_dir, ...
                                                       'constellation', constellation, 'num_looks', num_looks}];
                    end
                    switch num_spl_method
                        case 'auto'
                        OptionalArgs = [OptionalArgs, {'spline_method', spline_method}];
                        case 'manual'
                        OptionalArgs = [OptionalArgs, {'num_spl_manual', num_spl_manual}];
                    end
                    switch lambda_method
                        case 'manual'
                        OptionalArgs = [OptionalArgs, {'lambda_manual', lambda_manual}];
                    end
                    switch coll_proc
                        case 'prediction'
                        OptionalArgs = [OptionalArgs, {'coll_step_est', coll_step_est}];
                    end
                    if ismember(procType, {'temporal', 'temporal&NNI'})
                        modelConfig = phase_model_beta.loadConfig(app.ConfigRoot);
                        thresholdOptions = {
                            'min_period_days_method', 'min_period_days'
                            'min_coll_snr_method', 'min_coll_snr'
                            'min_coll_corr_samples_method', 'min_coll_corr_samples'
                            'spline_min_knot_intervals_method', 'spline_min_knot_intervals'
                            'spline_max_fraction_method', 'spline_max_fraction'
                        };
                        for thresholdIndex = 1:size(thresholdOptions,1)
                            methodName = thresholdOptions{thresholdIndex,1};
                            valueName = thresholdOptions{thresholdIndex,2};
                            if strcmp(modelConfig.(methodName), 'manual')
                                OptionalArgs = [OptionalArgs, ...
                                    {valueName, modelConfig.(valueName)}]; %#ok<AGROW>
                            end
                        end
                    end

                % - 5.2) Perform the time series modelling
                switch procType
                    case 'temporal'
                        [obs_p1, obs_p2, obs_p3, obs_p4, obs_p5] = ModellingInTime(...
                        detectedOS, outputDir, figsDir, dataIN_AOI, displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, xyIN_AOI, ...
                        varNoise_method, num_spl_method, lambda_method, coll_proc, OptionalArgs);

                        % convert variables names for report genartion
                        t_relTS = obs_p5{1, 1};
                        t_dateTS = obs_p5{1, 2};
                        displAOI_TS = zeros(size(obs_p5, 1), length(t_dateTS));
                        for i = 1:size(obs_p5, 1)
                            displAOI_TS(i, :) = obs_p5{i, 6}(1, :);
                        end

                        displAOI_stdTS = zeros(size(obs_p5, 1), length(t_dateTS));
                        for i = 1:size(obs_p5, 1)
                            displAOI_stdTS(i, :) = obs_p5{i, 6}(2, :);
                        end

                        % Update progress (70%)
                        app.notifyBetaProgress(70, 'Generating report variables...');
                        drawnow;

                    case 'temporal&NNI'
                        [obs_p1, obs_p2, obs_p3, obs_p4, obs_p5] = ModellingInTime(...
                        detectedOS, outputDir, figsDir, dataIN_AOI, displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, xyIN_AOI, ...
                        varNoise_method, num_spl_method, lambda_method, coll_proc, OptionalArgs);

                        % convert variables names for report genartion
                        t_relTS = obs_p5{1, 1};
                        t_dateTS = obs_p5{1, 2};
                        displAOI_TS = zeros(size(obs_p5, 1), length(t_dateTS));
                        for i = 1:size(obs_p5, 1)
                            displAOI_TS(i, :) = obs_p5{i, 6}(1, :);
                        end

                        displAOI_stdTS = zeros(size(obs_p5, 1), length(t_dateTS));
                        for i = 1:size(obs_p5, 1)
                            displAOI_stdTS(i, :) = obs_p5{i, 6}(2, :);
                        end

                        % perform spatial interpolation all over the AOI for each epoch
                        [displAOI_TS_NNI, displAOI_stdNNI, lonlatIN_AOI_NNI] = NaturalNeighborInterpolation(...
                            xyIN_AOI, displAOI_TS, displAOI_stdTS, xy_grid, centerline_data, figsDir, projDim, ...
                            t_dateTS, lonlat_grid_AOI, utmZone, markerSize);

                        % Update progress (70%)
                        app.notifyBetaProgress(70, 'Generating report variables...');
                        drawnow;
                end



                %% --- 6. Spatio-temporal displacement modelling ---
                % The displacement time series are processed all together to construct a
                % model that can be estimated at any position and time. This is chosen when
                % the goal is to have a continous and smooth surface that captures the
                % strongest deformation trends

                % - 6.1) Import the query coordiantes file
                xy_EXTR = [];
                xEXTR = [];
                yEXTR = [];
                latEXTR = [];
                lonEXTR = [];
                idEXTR = {};

                if flag_tsExtr

                    % Import the coordinates file
                    EXTR_table = readtable(filepath_EXTR, 'ReadVariableNames', true);

                    % Automatically detect the columns for coordinates
                    coord_columns = {};
                    for col = 1:width(EXTR_table)
                        % check if the column contains numeric data (latitude, longitude, or x/y)
                        if isnumeric(EXTR_table{:, col}) || islogical(EXTR_table{:, col})
                            % check column names for required coordinate names (lat, lon, x, y)
                            col_name = lower(EXTR_table.Properties.VariableNames{col});
                            if contains(col_name, 'lat') || contains(col_name, 'lon') || contains(col_name, 'x') || contains(col_name, 'y')
                                coord_columns{end+1} = EXTR_table.Properties.VariableNames{col};
                            end
                        end
                    end

                    if numel(coord_columns) == 2
                        % assuming the coordinates are in two columns
                        coords = EXTR_table{:, coord_columns};

                        % check if the columns are lat/lon or x/y and assign to appropriate variables
                        lowerColumns = lower(string(coord_columns));
                        latColumn = find(contains(lowerColumns,'lat'),1);
                        lonColumn = find(contains(lowerColumns,'lon'),1);
                        xColumn = find(lowerColumns=="x" | startsWith(lowerColumns,"x_"),1);
                        yColumn = find(lowerColumns=="y" | startsWith(lowerColumns,"y_"),1);
                        if ~isempty(latColumn) && ~isempty(lonColumn)
                            latEXTR = coords(:,latColumn);
                            lonEXTR = coords(:,lonColumn);
                        elseif ~isempty(xColumn) && ~isempty(yColumn)
                            % X/Y are explicitly in the selected project CRS.
                            xEXTR = coords(:,xColumn);
                            yEXTR = coords(:,yColumn);
                        else
                            error('Coordinates are not in lat/lon or x/y format.');
                        end
                    else
                        error('The table must contain exactly two coordinate columns (e.g., lat/lon or x/y).');
                    end

                    % check if the table contains the 'id' column, if not, set an artificial 'id' column
                    if ismember('id', EXTR_table.Properties.VariableNames)
                        idEXTR = EXTR_table.id;
                    else
                        idEXTR = cellstr(['q' num2str((1:height(EXTR_table))')]);
                    end

                    % check that EXTR coordinates are inside the AOI
                    if exist('latEXTR', 'var') && exist('lonEXTR', 'var')
                        [xEXTR, yEXTR] = phase_model_beta.projectWgs84( ...
                            projectCRS,[lonEXTR,latEXTR]);
                        inside_AOI_flag = inpolygon(xEXTR, yEXTR, xyAOI(:,1), xyAOI(:,2));
                    elseif exist('xEXTR', 'var') && exist('yEXTR', 'var')
                        inside_AOI_flag = inpolygon(xEXTR, yEXTR, xyAOI(:,1), xyAOI(:,2));
                        [latEXTR,lonEXTR] = projinv(projectCRS,xEXTR,yEXTR);
                    else
                        error('Coordinates not found or incorrect format.');
                    end

                    % check if there are coordinates outside the AOI
                    % all outside - error
                    if all(~inside_AOI_flag)
                        error('All coordinates are outside the AOI.');
                    end

                    % some outside - remove them
                    if any(~inside_AOI_flag)
                        warning('Some coordinates are outside the AOI. Removing them.');
                        latEXTR = latEXTR(inside_AOI_flag);
                        lonEXTR = lonEXTR(inside_AOI_flag);
                        xEXTR = xEXTR(inside_AOI_flag);
                        yEXTR = yEXTR(inside_AOI_flag);
                        if exist('idEXTR', 'var')
                            idEXTR = idEXTR(inside_AOI_flag);
                        end
                    end

                    % make the final vector
                    xy_EXTR = [xEXTR, yEXTR];

                end


                % Custom scripts are isolated from every automatic function.
                if isprop(app,'deterministicModelMode') && strcmp(app.deterministicModelMode,'custom')
                    customConfig = phase_model_beta.defaultConfig();
                    customNames = fieldnames(customConfig);
                    for customIndex = 1:numel(customNames)
                        customName = customNames{customIndex};
                        if isprop(app,customName), customConfig.(customName) = app.(customName); end
                    end
                    customConfig.stopCheck = @() phase_model_beta.throwIfStopped(app);
                % - 6.2) Perform the time series modelling
                switch procType
                    case 'spatialDET'
                        switch projDim
                            case '1D'
                                [lonlatIN_AOI_DET1D, t_dateTS, t_relTS, displAOI_TS, displAOI_TS_DET1D, stdAOI_DET1D, displEXTR] = ...
                                    STmodel_DET1D_custom(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, centerline_data, ...
                                    figsDir, minMonths, gS_input_path, gS_output_path, gS_job_path, gS_synth_path, ...
                                    detectedOS, varNoise_DET1D, varNoise_manual_DET1D, num_spl_method_DET1D, spline_method_DET1D, ...
                                    num_spl_row_manual_DET1D, num_spl_col_manual_DET1D, lambda_method_DET1D, lambda_manual_DET1D, ...
                                    utmZone, xyIN_AOI, step_t_ST, detrendMethodST, useInclinedMeansST, polyDegreeST, markerSize, xy_EXTR, customConfig, outputDir);
                            case '2D'
                                [lonlatIN_AOI_DET2D, t_dateTS, t_relTS, displAOI_TS, displAOI_TS_DET2D, stdAOI_DET2D, displEXTR] = ...
                                    STmodel_DET2D_custom(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, x_grid, y_grid, figsDir, minMonths, ...
                                    gS_input_path, gS_output_path, gS_job_path, detectedOS, varNoise_DET2D, varNoise_manual_DET2D, ...
                                    num_spl_method_DET2D, spline_method_DET2D, num_spl_row_manual_DET2D, num_spl_col_manual_DET2D, ...
                                    num_spl_t_manual_DET2D, lambda_method_DET2D, lambda_manual_DET2D, utmZone, xyIN_AOI, xyAOI, ...
                                    step_t_ST, detrendMethodST, useInclinedMeansST, polyDegreeST, markerSize, xy_EXTR, customConfig, outputDir);
                        end

                        % Update progress (70%)
                        app.notifyBetaProgress(70, 'Generating report variables...');
                        drawnow;

                    case 'spatialSTC'
                        switch projDim
                            case '1D'
                                [lonlatIN_AOI_STC1D, t_dateTS, t_relTS, displAOI_TS, displAOI_TS_STC1D, stdAOI_STC1D, displEXTR] = ...
                                    STmodel_STC1D_custom(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, centerline_data, figsDir, minMonths, ...
                                    gS_input_path, gS_output_path, gS_job_path, detectedOS, dtCov_STC1D, dsCov_STC1D, utmZone, xyIN_AOI, step_t_ST, ...
                                    detrendMethodST, useInclinedMeansST, polyDegreeST, tCovModel_STC2D, sCovModel_STC2D, markerSize, xy_EXTR, customConfig, outputDir);
                            case '2D'
                                [lonlatIN_AOI_STC2D, t_dateTS, t_relTS, displAOI_TS, displAOI_TS_STC2D, stdAOI_STC2D, displEXTR] = ...
                                    STmodel_STC2D_custom(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, x_grid, y_grid, useInclinedMeansST, figsDir, minMonths, ...
                                    gS_input_path, gS_output_path, gS_job_path, detectedOS, dtCov_STC2D, dsCov_STC2D, utmZone, xyIN_AOI, xyAOI, step_t_ST, ...
                                    detrendMethodST, polyDegreeST, tCovModel_STC2D, sCovModel_STC2D, markerSize, xy_EXTR, customConfig, outputDir);
                        end

                        % Update progress (70%)
                        app.notifyBetaProgress(70, 'Generating report variables...');
                        drawnow;

                end



                else
                % - 6.2) Perform the time series modelling
                switch procType
                    case 'spatialDET'
                        switch projDim
                            case '1D'
                                [lonlatIN_AOI_DET1D, t_dateTS, t_relTS, displAOI_TS, displAOI_TS_DET1D, stdAOI_DET1D, displEXTR] = ...
                                    STmodel_DET1D(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, centerline_data, ...
                                    figsDir, minMonths, gS_input_path, gS_output_path, gS_job_path, gS_synth_path, ...
                                    detectedOS, varNoise_DET1D, varNoise_manual_DET1D, num_spl_method_DET1D, spline_method_DET1D, ...
                                    num_spl_row_manual_DET1D, num_spl_col_manual_DET1D, lambda_method_DET1D, lambda_manual_DET1D, ...
                                    utmZone, xyIN_AOI, step_t_ST, detrendMethodST, useInclinedMeansST, polyDegreeST, markerSize, xy_EXTR);
                            case '2D'
                                [lonlatIN_AOI_DET2D, t_dateTS, t_relTS, displAOI_TS, displAOI_TS_DET2D, stdAOI_DET2D, displEXTR] = ...
                                    STmodel_DET2D(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, x_grid, y_grid, figsDir, minMonths, ...
                                    gS_input_path, gS_output_path, gS_job_path, detectedOS, varNoise_DET2D, varNoise_manual_DET2D, ...
                                    num_spl_method_DET2D, spline_method_DET2D, num_spl_row_manual_DET2D, num_spl_col_manual_DET2D, ...
                                    num_spl_t_manual_DET2D, lambda_method_DET2D, lambda_manual_DET2D, utmZone, xyIN_AOI, xyAOI, ...
                                    step_t_ST, detrendMethodST, useInclinedMeansST, polyDegreeST, markerSize, xy_EXTR);
                        end

                        % Update progress (70%)
                        app.notifyBetaProgress(70, 'Generating report variables...');
                        drawnow;

                    case 'spatialSTC'
                        switch projDim
                            case '1D'
                                [lonlatIN_AOI_STC1D, t_dateTS, t_relTS, displAOI_TS, displAOI_TS_STC1D, stdAOI_STC1D, displEXTR] = ...
                                    STmodel_STC1D(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, centerline_data, figsDir, minMonths, ...
                                    gS_input_path, gS_output_path, gS_job_path, detectedOS, dtCov_STC1D, dsCov_STC1D, utmZone, xyIN_AOI, step_t_ST, ...
                                    detrendMethodST, useInclinedMeansST, polyDegreeST, tCovModel_STC2D, sCovModel_STC2D, markerSize, xy_EXTR);
                            case '2D'
                                [lonlatIN_AOI_STC2D, t_dateTS, t_relTS, displAOI_TS, displAOI_TS_STC2D, stdAOI_STC2D, displEXTR] = ...
                                    STmodel_STC2D(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, x_grid, y_grid, useInclinedMeansST, figsDir, minMonths, ...
                                    gS_input_path, gS_output_path, gS_job_path, detectedOS, dtCov_STC2D, dsCov_STC2D, utmZone, xyIN_AOI, xyAOI, step_t_ST, ...
                                    detrendMethodST, polyDegreeST, tCovModel_STC2D, sCovModel_STC2D, markerSize, xy_EXTR);
                        end

                        % Update progress (70%)
                        app.notifyBetaProgress(70, 'Generating report variables...');
                        drawnow;

                end



                end

                %% --- 7. Variables for report generation ---
                % Prepare all the needed variables for the Excel report sheets

                % - 7.1) Define the variables of the final model
                switch procType
                    case 'temporal'
                        switch projDim
                            case '1D'
                                xyIN_AOI_FIN = xyIN_AOI;
                                lonlatIN_AOI_FIN = lonlatIN_AOI;
                                displAOI_FIN = displAOI_TS;
                                stdAOI_FIN = displAOI_stdTS;
                                displAOI_RAW = displIN_AOI - displAOI_TS(:,1);
                                t_relFIN = t_relTS; t_dateFIN = t_dateTS;
                                PSidIN_FIN = PSidIN_AOI;
                            case '2D'
                                xyIN_AOI_FIN = xyIN_AOI;
                                lonlatIN_AOI_FIN = lonlatIN_AOI;
                                displAOI_FIN = displAOI_TS;
                                stdAOI_FIN = displAOI_stdTS;
                                displAOI_RAW = displIN_AOI - displAOI_TS(:,1);
                                t_relFIN = t_relTS; t_dateFIN = t_dateTS;
                                PSidIN_FIN = PSidIN_AOI;
                        end
                    case 'temporal&NNI'
                        switch projDim
                            case '1D'
                                xyIN_AOI_FIN = centerline_data.xy_centerline;
                                lonlatIN_AOI_FIN = lonlatIN_AOI_NNI;
                                displAOI_FIN = displAOI_TS_NNI;
                                stdAOI_FIN = displAOI_stdNNI;
                                displAOI_RAW = displIN_AOI - displAOI_TS(:,1);
                                t_relFIN = t_relTS; t_dateFIN = t_dateTS;
                                PSidIN_FIN = (1:size(displAOI_FIN,1))';
                            case '2D'
                                xyIN_AOI_FIN = xy_grid;
                                lonlatIN_AOI_FIN = lonlatIN_AOI_NNI;
                                displAOI_FIN = displAOI_TS_NNI;
                                stdAOI_FIN = displAOI_stdNNI;
                                displAOI_RAW = displIN_AOI - displAOI_TS(:,1);
                                t_relFIN = t_relTS; t_dateFIN = t_dateTS;
                                PSidIN_FIN = (1:size(displAOI_FIN,1))';
                        end
                    case 'spatialDET'
                        switch projDim
                            case '1D'
                                xyIN_AOI_FIN = centerline_data.xy_centerline;
                                lonlatIN_AOI_FIN = lonlatIN_AOI_DET1D;
                                displAOI_FIN = displAOI_TS_DET1D;
                                stdAOI_FIN = stdAOI_DET1D;
                                displAOI_RAW = displIN_AOI - displAOI_TS(:,1);
                                t_relFIN = t_relTS; t_dateFIN = t_dateTS;
                                PSidIN_FIN = (1:size(displAOI_FIN,1))';
                            case '2D'
                                xyIN_AOI_FIN = xy_grid;
                                lonlatIN_AOI_FIN = lonlatIN_AOI_DET2D;
                                displAOI_FIN = displAOI_TS_DET2D;
                                stdAOI_FIN = stdAOI_DET2D;
                                displAOI_RAW = displIN_AOI - displAOI_TS(:,1);
                                t_relFIN = t_relTS; t_dateFIN = t_dateTS;
                                PSidIN_FIN = (1:size(displAOI_FIN,1))';
                        end
                    case 'spatialSTC'
                        switch projDim
                            case '1D'
                                xyIN_AOI_FIN = centerline_data.xy_centerline;
                                lonlatIN_AOI_FIN = lonlatIN_AOI_STC1D;
                                displAOI_FIN = displAOI_TS_STC1D;
                                stdAOI_FIN = stdAOI_STC1D;
                                displAOI_RAW = displIN_AOI - displAOI_TS(:,1);
                                t_relFIN = t_relTS; t_dateFIN = t_dateTS;
                                PSidIN_FIN = (1:size(displAOI_FIN,1))';
                            case '2D'
                                xyIN_AOI_FIN = xy_grid;
                                lonlatIN_AOI_FIN = lonlatIN_AOI_STC2D;
                                displAOI_FIN = displAOI_TS_STC2D;
                                stdAOI_FIN = stdAOI_STC2D;
                                displAOI_RAW = displIN_AOI - displAOI_TS(:,1);
                                t_relFIN = t_relTS; t_dateFIN = t_dateTS;
                                PSidIN_FIN = (1:size(displAOI_FIN,1))';
                        end
                end

                % - 7.2) Construct the information table for the raw displacements
                % create the headers
                headers1_RAW = {'ID', 'Latitude', 'Longitude', 'X (UTM)', 'Y (UTM)'};
                headers2_RAW = cellstr(datestr(t_dateIN, 'dd-mm-yyyy'));
                headers_RAW = [headers1_RAW, headers2_RAW'];

                % make the table contents
                table_RAW = [ ...
                    num2cell(PSidIN_AOI), ...              % PSid
                    num2cell(lonlatIN_AOI(:, 2)), ...      % latitude
                    num2cell(lonlatIN_AOI(:, 1)), ...      % longitude
                    num2cell(xyIN_AOI(:, 1)), ...          % X (UTM)
                    num2cell(xyIN_AOI(:, 2)), ...          % Y (UTM)
                    num2cell(displAOI_RAW) ...             % displacement at each epoch
                ];

                % combine headers and data
                fullTable_RAW = [headers_RAW', table_RAW'];


                % - 7.3) Figure of average RAW AOI displacement
                displINavg_RAW = mean(displAOI_RAW, 1);

                % perform linear regression to fit a linear trend line
                p = polyfit(t_relIN, displINavg_RAW, 1);
                trendLine = polyval(p, t_relIN);

                % trend velocity in mm/year
                slope_per_day = p(1);
                velocity_mm_per_year = slope_per_day * 365.25;

                % plot
                f = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
                plot(t_dateIN, displINavg_RAW, '-ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
                hold on
                plot(t_dateIN, trendLine, '-r', 'LineWidth', 1.5) % Trend line
                hold off
                xlim([t_dateIN(1), t_dateIN(end)]);
                xlabel('Time')
                ylabel('LOS displacement [mm]')
                set(gca, 'FontSize', 15)
                title(sprintf('Average scene displacement - %i PS', length(PSidIN_AOI)), 'FontSize', 20)
                legend('Mean displ.', 'Linear trend', 'Location', 'northeast')
                grid on
                annotation_text = sprintf('Trend velocity: %.2f mm/year', velocity_mm_per_year);
                x_pos = t_dateIN(1) + days(10);
                y_pos = max(displINavg_RAW) - 0.01 * range(displINavg_RAW);
                text(x_pos, y_pos, annotation_text, 'FontSize', 15, 'Color', 'red', 'FontWeight', 'bold')
                fig3_filename = strcat(figsDir, filesep, 'AOI_meanDispl_RAW.png');
                phase_model_beta.exportFigure(f,fig3_filename);
                close(f)


                % - 7.4) Construct the information table for the modelled displacements
                % initialize fields for the table
                num_PS = length(PSidIN_FIN);               % number of PS
                avg_velocity = zeros(num_PS, 1);           % average velocity in mm/year

                % compute average velocity for each PS
                for i = 1:num_PS
                    p = polyfit(t_relFIN, displAOI_FIN(i, :), 1);
                    avg_velocity(i) = p(1) * 365.25;       % slope (convert to mm/year)
                end

                % combine all data into a table
                headers_FIN = [{'ID new', 'Longitude', 'Latitude', 'X (UTM)', 'Y (UTM)', 'Avg. velocity [mm/yr]'}, ...
                    cellstr(datestr(t_dateFIN))'];

                % create the main table
                table_FIN = [num2cell(PSidIN_FIN), ...
                                 num2cell(lonlatIN_AOI_FIN), ...
                                 num2cell(xyIN_AOI_FIN), ...
                                 num2cell(avg_velocity), ...
                                 num2cell(displAOI_FIN)];

                % add headers
                fullTable_FIN = [headers_FIN', table_FIN'];


                % - 7.5) Make the figure of the interpolated average velocity over the AOI

                    % plot
                    f = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
                    geobasemap satellite;
                    hold on;
                    geoplot(lonlatAOI(:,2), lonlatAOI(:,1), 'k-', 'LineWidth', 2);
                    geoscatter(lonlatIN_AOI_FIN(:,2), lonlatIN_AOI_FIN(:,1), markerSize, avg_velocity, 'filled');
                    c = colorbar;
                    c.Label.String = 'Average LOS velocity [mm/year]';
                    set(gca, 'FontSize', 15)
                    title('Average velocity over AOI', 'FontSize', 20);
                    fig4_filename = fullfile(figsDir, 'AOI_AvgVelocity.png');
                    phase_model_beta.exportFigure(f,fig4_filename);
                    close(f);


                % - 7.6) Time-average model prediction uncertainty [mm]
                avg_std = mean(stdAOI_FIN,2,'omitnan');
                headers_STD = [{'ID new', 'Longitude', 'Latitude', ...
                    'X (UTM)', 'Y (UTM)', 'Mean prediction std [mm]'}, ...
                    cellstr(datestr(t_dateFIN))'];
                table_stdFIN = [num2cell(PSidIN_FIN), ...
                    num2cell(lonlatIN_AOI_FIN), num2cell(xyIN_AOI_FIN), ...
                    num2cell(avg_std), num2cell(stdAOI_FIN)];
                fullTable_stdFIN = [headers_STD',table_stdFIN'];

                % - 7.7) Make the figure of the uncertainties over the AOI

                    % plot
                    f = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
                    geobasemap satellite;
                    hold on;
                    geoplot(lonlatAOI(:,2), lonlatAOI(:,1), 'k-', 'LineWidth', 2);
                    geoscatter(lonlatIN_AOI_FIN(:,2), lonlatIN_AOI_FIN(:,1), markerSize, avg_std, 'filled');
                    c = colorbar;
                    c.Label.String = 'Mean model prediction std [mm]';
                    set(gca, 'FontSize', 15)
                    title('Time-average model prediction uncertainty over AOI', 'FontSize', 20);
                    fig5_filename = fullfile(figsDir, 'AOI_uncertainty.png');
                    phase_model_beta.exportFigure(f,fig5_filename);
                    close(f);


                % - 7.8) Descriptive anomalies of the full modelled field.
                % These scores are relative to the AOI, not hazard probabilities.
                anomaly = phase_model_beta.modelAnomalies( ...
                    displAOI_FIN,stdAOI_FIN);
                anomaly_table = table(PSidIN_FIN, ...
                    lonlatIN_AOI_FIN(:,1),lonlatIN_AOI_FIN(:,2), ...
                    xyIN_AOI_FIN(:,1),xyIN_AOI_FIN(:,2), ...
                    avg_velocity,anomaly.peakDeviationMm, ...
                    anomaly.peakContrast,100*anomaly.persistence, ...
                    string(anomaly.class), ...
                    'VariableNames', {'ID new','Longitude','Latitude', ...
                    'X (UTM)','Y (UTM)','Avg velocity [mm/yr]', ...
                    'Peak AOI deviation [mm]','Peak normalized contrast', ...
                    'Persistence [%]','Model anomaly'});

                % - 7.9) Shapefile creation and export
                % Check input consistency
                if size(lonlatIN_AOI_FIN, 1) ~= size(xyIN_AOI_FIN, 1) || ...
                   size(lonlatIN_AOI_FIN, 1) ~= size(PSidIN_FIN, 1) || ...
                   size(lonlatIN_AOI_FIN, 1) ~= size(avg_velocity, 1) || ...
                   size(lonlatIN_AOI_FIN, 1) ~= size(anomaly.peakContrast, 1) || ...
                   size(lonlatIN_AOI_FIN, 1) ~= size(displAOI_FIN, 1)
                    error('Input arrays have inconsistent sizes.');
                end

                % Initialize the structure array (vectorized)
                PHASEresults = struct('Geometry', repmat({'Point'}, size(lonlatIN_AOI_FIN, 1), 1), ...
                                   'X', num2cell(lonlatIN_AOI_FIN(:, 1)), ... % Longitude
                                   'Y', num2cell(lonlatIN_AOI_FIN(:, 2)), ... % Latitude
                                   'PSid', num2cell(PSidIN_FIN(:)), ...       % PS id
                                   'UTM_X', num2cell(xyIN_AOI_FIN(:, 1)), ... % UTM X
                                   'UTM_Y', num2cell(xyIN_AOI_FIN(:, 2)), ... % UTM Y
                                   'AvgVel', num2cell(avg_velocity), ...      % Average velocity
                                   'AvgStd', num2cell(avg_std), ...           % Average uncertainty
                                   'AnomPeak', num2cell(anomaly.peakContrast), ... % AOI-relative contrast
                                   'AnomClass', anomaly.class);         % Model anomaly

                % Create a FieldNameMap for .dbf field names
                FieldNameMap = containers.Map;
                FieldNameMap('PSid') = 'ID_new';
                FieldNameMap('UTM_X') = 'X_UTM';
                FieldNameMap('UTM_Y') = 'Y_UTM';
                FieldNameMap('AvgVel') = 'Avg_Vel';
                FieldNameMap('AvgStd') = 'Avg_Std';
                FieldNameMap('AnomPeak') = 'Anom_Peak';

                % Add displacement fields for each date
                for j = 1:size(displAOI_FIN, 2)
                    % Field name as Dyyyymmdd (e.g., 'D20190709')
                    fieldName = matlab.lang.makeValidName(['D' datestr(t_dateFIN(j), 'yyyymmdd')]);
                    % Assign displacement values directly from displAOI_FIN
                    v = displAOI_FIN(:, j);
                    for i = 1:length(PHASEresults)
                        PHASEresults(i).(fieldName) = v(i);
                    end
                end

                % Export to shapefile
                try
                    shapewrite(PHASEresults, filenameOUT_s);
                catch ME
                    error('Failed to write shapefile: %s', ME.message);
                end

                % Setting the CRS (WGS84, EPSG:4326)
                try
                    proj = geocrs(4326, 'Authority', 'EPSG');
                    strWkt = wktstring(proj, 'Version', 'wkt1');
                catch
                    % Fallback WKT for WGS84
                    strWkt = ['GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],' ...
                              'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'];
                end

                % Create the .prj file
                prjFile = fullfile(filesDir, 'shp', strcat(filenameOUT, '.prj'));
                fid = fopen(prjFile, 'w');
                if fid == -1
                    error('Failed to create .prj file: %s', prjFile);
                end
                fprintf(fid, '%s', strWkt);
                fclose(fid);


                % - 7.10) Export the .mat file with the same structure as shapefile
                save(fullfile(filesDir, 'mat', 'PHASEresults.mat'), 'PHASEresults');


                % - 7.11) Export multi-band GeoTIFF for 2D gridded models
                if strcmp(projDim, '2D') && ~strcmp(procType,'temporal')
                    fprintf('Exporting multi-band GeoTIFF of the 2D spatio-temporal model...\n');

                    % 7.11.1) Dimensions
                    n_rows = size(x_grid, 1);
                    n_cols = size(x_grid, 2);
                    n_times = size(displAOI_FIN, 2);

                    % 7.11.2) Pre-allocate 3D rasters (NaN automatically acts as NoData outside the AOI)
                    tiff_matrix_displ = NaN(n_rows, n_cols, n_times);
                    tiff_matrix_std   = NaN(n_rows, n_cols, n_times);

                    % 7.11.3) Fill the matrix epoch by epoch
                    for t = 1:n_times
                        % Displacement slice
                        temp_slice_displ = NaN(n_rows, n_cols);
                        temp_slice_displ(inAOI_xy) = displAOI_FIN(:, t);
                        tiff_matrix_displ(:,:,t) = flip(temp_slice_displ, 1);

                        % Uncertainty slice
                        temp_slice_std = NaN(n_rows, n_cols);
                        temp_slice_std(inAOI_xy) = stdAOI_FIN(:, t);
                        tiff_matrix_std(:,:,t) = flip(temp_slice_std, 1);
                    end

                    % 7.11.4) Create the spatial reference object
                    pixel_size_x = abs(x_grid(1,2) - x_grid(1,1));
                    pixel_size_y = abs(y_grid(2,1) - y_grid(1,1));

                    x_limits = [min(x_grid(:)) - pixel_size_x/2, max(x_grid(:)) + pixel_size_x/2];
                    y_limits = [min(y_grid(:)) - pixel_size_y/2, max(y_grid(:)) + pixel_size_y/2];

                    R = maprefcells(x_limits, y_limits, [n_rows, n_cols]);

                    % 7.11.5) Write rasters in the selected project CRS
                    filenameOUT_tiff_displ = fullfile(filesDir, strcat(filenameOUT, '_model.tif'));
                    filenameOUT_tiff_std   = fullfile(filesDir, strcat(filenameOUT, '_uncertainty.tif'));

                    geotiffwrite(filenameOUT_tiff_displ,tiff_matrix_displ,R, ...
                        'CoordRefSysCode',epsgCode);
                    geotiffwrite(filenameOUT_tiff_std,tiff_matrix_std,R, ...
                        'CoordRefSysCode',epsgCode);
                    fprintf('GeoTIFFs exported with EPSG:%d.\n',epsgCode);
                    end

                pause(10)
                app.notifyBetaProgress(80, 'Extrapolating time series...');
                drawnow;



                %% --- 8. Time series interpolation ---
                % Extrapolate the displacement time series in correspondence of the
                % provided coordinates.

                if flag_tsExtr && ~isempty(xy_EXTR)

                    % - 8.1) Bridge variables for interpolation
                    nEXTR_points = size(xy_EXTR, 1);
                    displEXTR_cell = cell(nEXTR_points, 1);
                    EXTR_flag = zeros(nEXTR_points, 1);

                    if contains(procType, 'spatial')
                        % Map native matrix output from STmodels to cell array
                        for i = 1:nEXTR_points
                            displEXTR_cell{i} = displEXTR(i, :);
                            EXTR_flag(i) = 1; % Mark as valid
                        end
                    else
                        % Fallback snapping for 'temporal' or 'temporal&NNI' models
                        for i = 1:nEXTR_points
                            % Use xy_EXTR directly to guarantee variable scope
                            distances = sqrt((xyIN_AOI_FIN(:, 1) - xy_EXTR(i, 1)).^2 + (xyIN_AOI_FIN(:, 2) - xy_EXTR(i, 2)).^2);
                            [min_distance, idx] = min(distances);

                            % Safe fallback for pure temporal models where resolutions might be 0 or NaN
                            if strcmp(projDim, '1D')
                                tolerance = max([cline_resolution, 50], [], 'omitnan');
                            else
                                tolerance = max([grid_resolution, 50], [], 'omitnan');
                            end

                            if min_distance > tolerance
                                warning('No PS found near query point %s within tolerance (%.1f m).', idEXTR{i}, min_distance);
                                EXTR_flag(i) = NaN;
                                displEXTR_cell{i} = NaN;
                            else
                                EXTR_flag(i) = idx;
                                displEXTR_cell{i} = displAOI_FIN(idx, :);
                            end
                        end
                    end

                    % Reassign back to displEXTR for 8.3 plotting and Step 10 Excel export
                    displEXTR = displEXTR_cell;


                    % - 8.2) Figures creation for each extrapolated time series
                    % initialize a cell array to store filenames
                    figEXTRnames = cell(nEXTR_points, 1);

                    for i = 1:nEXTR_points
                        if isnan(EXTR_flag(i))
                            continue;    % skip unmatched points
                        end

                        % plot
                        f = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
                        plot(t_dateFIN, displEXTR{i}, '-ko', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k');
                        grid on;
                        xlabel('Time');
                        ylabel('LOS displacement [mm]');
                        set(gca, 'FontSize', 15);
                        title(sprintf('Modelled displacement for point %s (%.3f, %.3f)', idEXTR{i}, latEXTR(i), lonEXTR(i)), 'Interpreter', 'none', 'FontSize', 20);
                        xtickangle(45);
                        figEXTRnames{i} = fullfile(figsDir, sprintf('ExtrDispl_%s.png', idEXTR{i}));
                        print(f, figEXTRnames{i}, '-dpng', '-r300');
                        close(f);
                    end

                end

                % Update progress (90%)
                app.notifyBetaProgress(90, 'Generating Excel report...');
                drawnow;



                %% --- 10. Excel report ---
                % Creation and export of processed data and analysis results in an Excel
                % spreadsheet for easy interpretability

                % - 10.3) Prepare and fill the first sheet: 'General'
                sheetName1 = 'General';

                % general variables
                OUTtitle = sprintf('%s PHASE geospatial model for %s', projDim, municipality);
                OUTsubtitle = 'Automatic report generated by PHASE';
                OUTmunicipality = sprintf('Municipality: %s', municipality);
                OUTcountry = sprintf('Country: %s', country);
                OUTcurrentDT = sprintf('%s', currentDT);
                OUTprojDim = projDim;
                OUTreflat = num2str(sprintf('%.8f', ref_centre_lonlat(2)));
                OUTreflon = num2str(sprintf('%.8f', ref_centre_lonlat(1)));
                OUTrefRadius = sprintf('%i m', ref_radius);
                OUTnumPS = num2str(size(displIN, 1));
                OUTnumPS_AOI = num2str(size(displIN_AOI, 1));
                OUTobsPeriod = sprintf('%s - %s', t_dateIN(1), t_dateIN(end));
                OUTdays = sprintf('%s days', num2str(t_relIN(end)));
                OUTstampsParamsAvailable = isfile(filepathREF);

                % create cell array for the left side of the sheet
                dataSheet1 = {
                    OUTtitle, '';                                      % title
                    OUTsubtitle, '';                                   % subtitle
                    ['Created on ', OUTcurrentDT], '';
                    '', '';                                            % blank row
                    OUTmunicipality, '';                               % municipality
                    OUTcountry, '';                                    % country
                    '', '';                                            % blank row
                    ['Chosen modelling: ', OUTprojDim], '';
                    ['Project CRS: EPSG:', num2str(epsgCode)], '';            % project spatial dimension
                    '', '';                                            % blank row
                    'Unwrapping parameters', '';
                    ['Reference latitude [deg]: ', OUTreflat], '';     % reference latitude
                    ['Reference longitude [deg]: ', OUTreflon], '';    % reference longitude
                    ['Reference Radius: ', OUTrefRadius], '';      % reference radius
                    '', '';                                            % blank row
                    'Persistent Scatterers', '';
                    ['Total PS: ', OUTnumPS], '';                      % total Persistent Scatterers
                    ['PS in AOI: ', OUTnumPS_AOI], '';                 % Persistent Scatterers in AOI
                    '', '';                                            % blank row
                    'Project time span', '';
                    ['Observation Period: ', OUTobsPeriod], '';        % observation period
                    ['Observation Length: ', OUTdays], '';             % length in days
                };

                % handle StaMPS parameters if available
                if OUTstampsParamsAvailable
                    % load StaMPS parameters
                    stampsParams = load(filepathREF);

                    % convert StaMPS parameters to a cell array of name-value pairs
                    paramsData = prepareStaMPSParams(stampsParams);

                    % avoid to expand large vectors
                    for i = 1:size(paramsData, 1)
                        val = paramsData{i, 2};

                        if (isnumeric(val) || islogical(val)) && numel(val) > 1
                            if numel(val) > 10
                                % if big write a summary
                                paramsData{i, 2} = sprintf('[Array %dx%d %s]', size(val,1), size(val,2), class(val));
                            else
                                % if small convert to text
                                paramsData{i, 2} = mat2str(val);
                            end
                        elseif ischar(val) && length(val) > 30000
                            % safety 1
                            paramsData{i, 2} = '[Too long text for Excel]';
                        elseif iscell(val) || isstruct(val)
                            % safety 2
                            paramsData{i, 2} = sprintf('[%s data]', class(val));
                        end
                    end

                    % add StaMPS parameters to the left-side data
                    dataSheet1 = [dataSheet1; {'', ''}; {'StaMPS processing parameters:', ''}; paramsData]; % Append parameters
                end

                % write data to the Excel file
                writecell(dataSheet1, filenameOUT_e, 'Sheet', sheetName1, 'Range', 'A1');

                % full command to call Python script
                pythonScript = fullfile('pythonScripts', 'reportSheet1.py');
                pythonExecutable = pyenv().Executable;

                % define figure paths and positions
                figures = {logo_filename, fig2_filename, fig1_filename};
                positions = {'F1', 'E8', 'E30'};
                dimensions = {[366, 180], [900, 450], [900, 450]};

                % Combine figure details into a single string
                figuresAndPositions = '';
                for i = 1:length(figures)
                    fig_details = sprintf('%s,%s,%d,%d', figures{i}, positions{i}, dimensions{i}(1), dimensions{i}(2));
                    if i < length(figures)
                        figuresAndPositions = strcat(figuresAndPositions, fig_details, ';');
                    else
                        figuresAndPositions = strcat(figuresAndPositions, fig_details);
                    end
                end

                % prepare the command
                command = sprintf('%s "%s" "%s" "%s" "%s"', pythonExecutable, pythonScript, filenameOUT_e, figuresAndPositions, sheetName1);

                % run the Python script from MATLAB
                [status, cmdout] = phase_model_beta.runCommandHidden(app, command, 'Excel report formatting');

                % check for errors
                if status == 0
                    disp('Python script for Excel sheet 1 executed successfully.');
                    disp(cmdout);
                else
                    disp('Error executing Python script for Excel sheet 1:');
                    disp(cmdout);
                end


                % - 10.4) Prepare and fill the second sheet: 'Raw displacement'
                sheetName2 = 'Raw displacement';

                % define the placeholders for the title, subtitles, and figure space
                title_s2 = 'Raw displacement time series';
                subtitle1_s2 = ['Each PS displacement time series is obtained subtracting the value ' ...
                    'of the modelled trend displacement at the first epoch, to the whole time series.'];
                subtitle2_s2 = 'Coordinates are WGS84 [deg] and project UTM [m]; see the EPSG code in General.';
                subtitle3_s2 = 'Displacement is given along the satellite Line Of Sight (LOS) [mm].';
                subtitle4_s2 = 'Persistent Scatterer data - inside the AOI';

                % combine all data into a cell array for the second sheet
                dataSheet2 = {
                    title_s2, '';            % title
                    subtitle1_s2, '';        % first subtitle
                    subtitle2_s2, '';        % second subtitle
                    subtitle3_s2, '';        % third subtitle
                    '', '';                  % blank row
                };

                % add blank rows for the figure space
                figureSpace = cell(25, 2);   % 25 empty rows for the figure
                dataSheet2 = [dataSheet2; figureSpace];

                % add the third subtitle
                dataSheet2 = [dataSheet2; {subtitle4_s2, ''}];

                % write the placeholders to the Excel file
                writecell(dataSheet2, filenameOUT_e, 'Sheet', sheetName2, 'Range', 'A1');

                % write the table below the placeholders
                startTableRow = size(dataSheet2, 1) + 1;
                startTableCell = sprintf('A%d', startTableRow);
                writecell(fullTable_RAW, filenameOUT_e, 'Sheet', sheetName2, 'Range', startTableCell);

                disp('Second sheet placeholders and table written successfully.');

                % define the Python script and arguments
                pythonScript = fullfile('pythonScripts', 'reportSheet2.py');
                pythonExecutable = pyenv().Executable;

                % command to format the second sheet and execution
                command = sprintf('%s "%s" "%s" "%s" "%s"', pythonExecutable, pythonScript, filenameOUT_e, fig3_filename, sheetName2);
                [status, cmdout] = phase_model_beta.runCommandHidden(app, command, 'Excel report formatting');

                % check for success
                if status == 0
                    disp('Python script for formatting the second sheet executed successfully.');
                    disp(cmdout);
                else
                    disp('Error executing Python script for the second sheet:');
                    disp(cmdout);
                end


                % - 10.5) Prepare and fill the third sheet: 'Modelled displacement'
                sheetName3 = 'Modelled displacement';

                % define the placeholders for the title, subtitles, and figure space
                title_s3 = 'Modelled displacement time series';
                subtitle1_s3 = ['Each PS displacement modelled time series is obtained through the ' ...
                    'contribution of a cubic splines deterministic trend and a stochastic component.'];
                if exist('customConfig','var')
                    subtitle1_s3 = ['Custom two-stage model: assigned/estimated deterministic contribution ' ...
                        'plus residual spline (DET) or collocation (STC). See custom_model.json for formula and parameters.'];
                end
                subtitle2_s3 = 'Coordinates are WGS84 [deg] and project UTM [m]; see the EPSG code in General.';
                subtitle3_s3 = 'Displacement is given along the satellite Line Of Sight (LOS) [mm]. The average velocity is in [mm/year].';
                subtitle4_s3 = 'Persistent Scatterer modelled data - inside the AOI';

                % combine all data into a cell array for the third sheet
                dataSheet3 = {
                    title_s3, '';            % title
                    subtitle1_s3, '';        % first subtitle
                    subtitle2_s3, '';        % second subtitle
                    subtitle3_s3, '';        % third subtitle
                    '', '';                  % blank row
                };

                % add blank rows for the figure space
                figureSpace = cell(25, 2);   % 25 empty rows for the figure
                dataSheet3 = [dataSheet3; figureSpace];

                % add the third subtitle
                dataSheet3 = [dataSheet3; {subtitle4_s3, ''}];

                % write the placeholders to the Excel file
                writecell(dataSheet3, filenameOUT_e, 'Sheet', sheetName3, 'Range', 'A1');

                % write the table below the placeholders
                startTableRow = size(dataSheet3, 1) + 1;
                startTableCell = sprintf('A%d', startTableRow);
                writecell(fullTable_FIN, filenameOUT_e, 'Sheet', sheetName3, 'Range', startTableCell);

                disp('Third sheet placeholders and table written successfully.');

                % define the Python script and arguments
                pythonScript = fullfile('pythonScripts', 'reportSheet3.py');
                pythonExecutable = pyenv().Executable;

                % command to format the third sheet and execution
                command = sprintf('%s "%s" "%s" "%s" "%s"', pythonExecutable, pythonScript, filenameOUT_e, fig4_filename, sheetName3);
                [status, cmdout] = phase_model_beta.runCommandHidden(app, command, 'Excel report formatting');

                % check for success
                if status == 0
                    disp('Python script for formatting the third sheet executed successfully.');
                    disp(cmdout);
                else
                    disp('Error executing Python script for the third sheet:');
                    disp(cmdout);
                end


                % - 10.6) Prepare and fill the fourth sheet: 'Uncertainty'
                sheetName4 = 'Uncertainty';

                % define the placeholders for the title, subtitles, and figure space
                title_s4 = 'Uncertainty of modelled displacement time series';
                if strcmp(procType,'temporal&NNI')
                    subtitle1_s4 = ['NNI interpolates the temporal model standard deviations ' ...
                        'at grid locations; no interpolation covariance is propagated.'];
                else
                    subtitle1_s4 = ['Per-epoch prediction standard deviation follows the final ' ...
                        'modelling stage; the summary is its arithmetic mean over time.'];
                end
                subtitle2_s4 = 'Coordinates are WGS84 [deg] and project UTM [m]; see the EPSG code in General.';
                if exist('customConfig','var')
                    subtitle1_s4 = ['Custom two-stage local uncertainty includes parameter/residual cross covariance and temporal referencing. ' ...
                        'Conditional on covariance/spline settings and outlier mask; excludes hyperparameter uncertainty, spline bias and model discrepancy.'];
                end
                subtitle3_s4 = 'Prediction standard deviation is in LOS millimetres [mm].';
                subtitle4_s4 = 'Persistent Scatterer modelled data uncertainty - inside the AOI';

                % combine all data into a cell array for the third sheet
                dataSheet4 = {
                    title_s4, '';            % title
                    subtitle1_s4, '';        % first subtitle
                    subtitle2_s4, '';        % second subtitle
                    subtitle3_s4, '';        % third subtitle
                    '', '';                  % blank row
                };

                % add blank rows for the figure space
                figureSpace = cell(25, 2);   % 25 empty rows for the figure
                dataSheet4 = [dataSheet4; figureSpace];

                % add the third subtitle
                dataSheet4 = [dataSheet4; {subtitle4_s4, ''}];

                % write the placeholders to the Excel file
                writecell(dataSheet4, filenameOUT_e, 'Sheet', sheetName4, 'Range', 'A1');

                % write the table below the placeholders
                startTableRow = size(dataSheet4, 1) + 1;
                startTableCell = sprintf('A%d', startTableRow);
                writecell(fullTable_stdFIN, filenameOUT_e, 'Sheet', sheetName4, 'Range', startTableCell);

                disp('Fourth sheet placeholders and table written successfully.');

                % define the Python script and arguments
                pythonScript = fullfile('pythonScripts', 'reportSheet4.py');
                pythonExecutable = pyenv().Executable;

                % command to format the third sheet and execution
                command = sprintf('%s "%s" "%s" "%s" "%s"', pythonExecutable, pythonScript, filenameOUT_e, fig5_filename, sheetName4);
                [status, cmdout] = phase_model_beta.runCommandHidden(app, command, 'Excel report formatting');

                % check for success
                if status == 0
                    disp('Python script for formatting the fourth sheet executed successfully.');
                    disp(cmdout);
                else
                    disp('Error executing Python script for the fourth sheet:');
                    disp(cmdout);
                end


                % - 10.7) Prepare and fill the fifth sheet: 'Alerts'
                sheetName5 = 'AOI anomalies';

                % define the placeholders for the title and subtitles
                title_s5 = 'Modelled deformation anomalies relative to the AOI';
                subtitle1_s5 = ['The full modelled time series is compared with the ' ...
                    'AOI median at each epoch. Scores are descriptive, not hazard probabilities.'];
                subtitle2_s5 = 'Persistent/transient labels require two adjacent epochs above three robust AOI scales.';
                subtitle3_s5 = 'Modelled positions inside the AOI';

                % combine all data into a cell array for the fourth sheet
                dataSheet5 = {
                    title_s5, '';            % title
                    subtitle1_s5, '';        % first subtitle
                    subtitle2_s5, '';        % second subtitle
                    '', '';                  % blank row
                    subtitle3_s5, '';        % third subtitle
                };

                % write the placeholders to the Excel file
                writecell(dataSheet5, filenameOUT_e, 'Sheet', sheetName5, 'Range', 'A1');

                % write the table below the titles
                startTableRow = size(dataSheet5, 1) + 1;
                startTableCell = sprintf('A%d', startTableRow);
                writetable(anomaly_table, filenameOUT_e, 'Sheet', sheetName5, 'Range', startTableCell);

                disp('Fifth sheet placeholders and table written successfully.');

                % define the Python script and arguments
                pythonScript = fullfile('pythonScripts', 'reportSheet5.py');
                pythonExecutable = pyenv().Executable;

                % command to format the fourth sheet and execution
                command = sprintf('%s "%s" "%s" "%s"', pythonExecutable, pythonScript, filenameOUT_e, sheetName5);
                [status, cmdout] = phase_model_beta.runCommandHidden(app, command, 'Excel report formatting');

                % check for success
                if status == 0
                    disp('Python script for formatting the fourth sheet executed successfully.');
                    disp(cmdout);
                else
                    disp('Error executing Python script for the fifth sheet:');
                    disp(cmdout);
                end


                % - 10.8) Prepare and fill the fifth sheet: 'Interpolation'
                if flag_tsExtr && ~all(isnan(EXTR_flag(:)))

                    sheetName6 = 'Interpolation';

                    % define the placeholders for the title and subtitles
                    title_s6 = 'Time series interpolation';
                    if contains(procType, 'spatial')
                        subtitle1_s6 = 'Displacement time series are extrapolated by rigorously evaluating the continuous spatio-temporal model at the exact input coordinates.';
                    else
                        subtitle1_s6 = 'Displacement time series are estimated by snapping to the nearest valid observation or grid node within a safe tolerance radius.';
                    end
                    subtitle2_s6 = 'IDs are preserved when given and automatically generated otherwise.';
                    subtitle3_s6 = 'Interpolated time series - data';

                    % combine all data into a cell array for the fifth sheet
                    dataSheet6 = {
                        title_s6, '';            % title
                        subtitle1_s6, '';        % first subtitle
                        subtitle2_s6, '';        % second subtitle
                        '', '';                  % blank row
                        subtitle3_s6, '';        % third subtitle
                    };

                    % write the placeholders to the Excel file
                    writecell(dataSheet6, filenameOUT_e, 'Sheet', sheetName6, 'Range', 'A1');

                    % define the starting row for writing time series data
                    startRow = size(dataSheet6, 1) + 1;

                    % loop through each extracted point
                    for i = 1:nEXTR_points
                        if isnan(EXTR_flag(i))
                            continue;    % skip unmatched points
                        end

                        % define headers for the time series data
                        headers = [{'ID', 'Latitude', 'Longitude', 'X (UTM)', 'Y (UTM)'}, cellstr(datestr(t_dateFIN))'];

                        % retrieve coordinate information
                        coordRow = {idEXTR{i}, latEXTR(i), lonEXTR(i), xEXTR(i), yEXTR(i)};

                        % merge headers, coordinate row, and time series data
                        fullData = [headers; coordRow, num2cell(displEXTR{i})];

                        % write to Excel
                        rangee = sprintf('A%d', startRow);
                        writecell(fullData, filenameOUT_e, 'Sheet', sheetName6, 'Range', rangee);

                        % leave 25 rows empty for the figure
                        startRow = startRow + size(fullData, 1) + 25;
                    end

                    disp('Sixth sheet placeholders and table written successfully.');

                    % define the Python script and arguments
                    pythonScript = fullfile('pythonScripts', 'reportSheet6.py');
                    pythonExecutable = pyenv().Executable;

                    % convert figFilenames to a comma-separated string
                    figuresArg = strjoin(figEXTRnames, ',');

                    % command to format the sixth sheet and execution
                    command = sprintf('%s "%s" "%s" "%s" "%s"', pythonExecutable, pythonScript, filenameOUT_e, sheetName6, figuresArg);
                    [status, cmdout] = phase_model_beta.runCommandHidden(app, command, 'Excel report formatting');

                    % check for success
                    if status == 0
                        disp('Python script for formatting the sixth sheet executed successfully.');
                        disp(cmdout);
                    else
                        disp('Error executing Python script for the sixth sheet:');
                        disp(cmdout);
                    end

                end


                % - 10.9) Prepare and fill the sheet: 'Modelled PS' (Optional)
                % Checking if the user checked the UI box
                if flag_PSinterp

                    sheetName7 = 'Modelled PS';

                    % The modelled displacement at the exact observation points is stored in displAOI_TS
                    title_s7 = 'Modelled displacement at observation points';
                    subtitle1_s7 = 'This sheet contains the evaluation of the model exactly at the input PS coordinates.';

                    dataSheet7 = {
                        title_s7, '';
                        subtitle1_s7, '';
                        '', '';
                        };

                    % Write the title placeholders
                    writecell(dataSheet7, filenameOUT_e, 'Sheet', sheetName7, 'Range', 'A1');

                    % --- THE FIX: Select the correct time vector to match displAOI_TS ---
                    if contains(procType, 'spatial')
                        obs_dates = t_dateIN;
                    else
                        obs_dates = t_dateTS;
                    end

                    % Build the data table matching your RAW format
                    headers_OBS = [{'ID', 'Latitude', 'Longitude', 'X (UTM)', 'Y (UTM)'}, cellstr(datestr(obs_dates))'];

                    table_OBS = [num2cell(PSidIN_AOI), ...
                        num2cell(lonlatIN_AOI(:, 2)), ...           % Latitude
                        num2cell(lonlatIN_AOI(:, 1)), ...           % Longitude
                        num2cell(xyIN_AOI(:, 1)), ...               % X (UTM)
                        num2cell(xyIN_AOI(:, 2)), ...               % Y (UTM)
                        num2cell(displAOI_TS - displAOI_TS(:,1))];  % Modeled series

                    fullTable_OBS = [headers_OBS', table_OBS'];

                    % Write the table below the titles
                    startTableRow = size(dataSheet7, 1) + 1;
                    writecell(fullTable_OBS, filenameOUT_e, 'Sheet', sheetName7, 'Range', sprintf('A%d', startTableRow));

                    % define the Python script and arguments
                    pythonScript = fullfile('pythonScripts', 'reportSheet7.py');
                    pythonExecutable = pyenv().Executable;

                    % command to format the sixth sheet and execution
                    command = sprintf('%s "%s" "%s" "%s"', pythonExecutable, pythonScript, filenameOUT_e, sheetName7);
                    [status, cmdout] = phase_model_beta.runCommandHidden(app, command, 'Excel report formatting');

                    % check for success
                    if status == 0
                        disp('Python script for formatting the seventh sheet executed successfully.');
                        disp(cmdout);
                    else
                        disp('Error executing Python script for the seventh sheet:');
                        disp(cmdout);
                    end

                    disp('Modeled PS exported to Excel successfully.');
                end

                % update progress (100%)
                app.notifyBetaProgress(100, 'Processing complete!');
                drawnow;


                fprintf('PHASE: end of processing. \n')


                % Progress is shown in the PHASE Model run monitor.

            catch err

                % 1. Extract the exact line and function where the error occurred
                if ~isempty(err.stack)
                    errorLine = err.stack(1).line;
                    errorFunc = err.stack(1).name;
                    % Build a detailed message
                    detailedMessage = sprintf('Error in function: %s (Line %d)\n\n%s', ...
                        errorFunc, errorLine, err.message);
                else
                    detailedMessage = err.message;
                end

                % 3. update UI label

                % 4. hand the error to the visible standalone controller
                app.notifyBetaProgress(100, ['Processing failed: ' err.message]);
                rethrow(err);

            end

        end
    end
end
