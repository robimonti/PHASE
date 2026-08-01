"""Extract PHASE_model.mlapp into a reproducible editable MATLAB class."""

from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[1]
MLAPP = ROOT / "legacy" / "PHASE_model.mlapp"
OUTPUT = ROOT / "+phase_model_beta" / "LegacyEngine.m"


def extract(xml: str) -> str:
    cdata_start = xml.index("<![CDATA[") + len("<![CDATA[")
    cdata_end = xml.rindex("]]>")
    code = xml[cdata_start:cdata_end].replace("\r\n", "\n")
    code = code.replace(
        "classdef PHASE_model < matlab.apps.AppBase",
        "classdef LegacyEngine < matlab.apps.AppBase",
        1,
    )
    code = code.replace(
        "function app = PHASE_model",
        "function app = LegacyEngine",
        1,
    )
    code = code.replace("properties (Access = private)", "properties (Access = public)")
    code = code.replace("methods (Access = private)", "methods (Access = public)")
    code = code.replace(
        "fileparts(mfilename('fullpath'))",
        "phase_model_beta.projectRoot()",
    )
    code = code.replace(
        'fileparts(mfilename("fullpath"))',
        "phase_model_beta.projectRoot()",
    )
    startup_anchor = "            addpath(fullfile(pwd, 'MatlabFunctions'));"
    if code.count(startup_anchor) < 1:
        raise RuntimeError("Could not locate PHASE Model startup path anchor")
    code = code.replace(
        startup_anchor,
        "            rootDir = phase_model_beta.projectRoot();\n"
        "            cd(rootDir);\n"
        "            addpath(fullfile(rootDir, 'MatlabFunctions'));\n"
        "            phase_model_beta.themeLegacyEngine(app);\n"
        "            app.pythonPath = phase_model_beta.resolvePythonPath(app.pythonPath);\n"
        "            app.pythoninstallationpathEditField.Value = app.pythonPath;",
        1,
    )
    code = code.replace(
        startup_anchor,
        "                addpath(fullfile(phase_model_beta.projectRoot(), 'MatlabFunctions'));",
    )
    code = code.replace(
        "fullfile(pwd, 'PHASE_logo.png')",
        "fullfile(phase_model_beta.projectRoot(), 'PHASE_logo.png')",
    )

    run_root_anchor = "                % --- 0. Prepare the environment ---"
    if code.count(run_root_anchor) != 1:
        raise RuntimeError("Could not locate Model processing-root anchor")
    code = code.replace(
        run_root_anchor,
        """\
                % Every legacy relative path is rooted explicitly for the
                % complete run. Results are stored beside the visible PHASE
                % shortcuts rather than inside the editable engine clone.
                runtimeRoot = phase_model_beta.projectRoot();
                previousRunFolder = pwd;
                runFolderCleanup = onCleanup(@() cd(previousRunFolder)); %#ok<NASGU>
                cd(runtimeRoot);
                outputRoot = fileparts(runtimeRoot);

                % --- 0. Prepare the environment ---
""",
        1,
    )

    output_scan_anchor = """\
                allFolders = dir();
                allFolderNames = {allFolders([allFolders.isdir]).name};
                outputFolders = allFolderNames(startsWith(allFolderNames, baseFolderName));
"""
    if code.count(output_scan_anchor) != 2:
        raise RuntimeError("Could not locate both Model output-folder scans")
    code = code.replace(
        output_scan_anchor,
        """\
                allFolders = dir(outputRoot);
                allFolderNames = {allFolders([allFolders.isdir]).name};
                outputFolders = allFolderNames(startsWith(allFolderNames, baseFolderName));
""",
    )

    output_remove_anchor = """\
                            folderToRemove = outputFolders{i};
                            if ~strcmp(folderToRemove, '.') && ~strcmp(folderToRemove, '..')
                                rmdir(folderToRemove, 's');
"""
    if code.count(output_remove_anchor) != 1:
        raise RuntimeError("Could not locate Model output-folder removal")
    code = code.replace(
        output_remove_anchor,
        """\
                            folderToRemove = outputFolders{i};
                            if ~strcmp(folderToRemove, '.') && ~strcmp(folderToRemove, '..')
                                rmdir(fullfile(outputRoot,folderToRemove), 's');
""",
        1,
    )

    output_create_anchor = """\
                % create the folder
                mkdir(outputDir);"""
    if code.count(output_create_anchor) != 1:
        raise RuntimeError("Could not locate Model output-folder creation")
    code = code.replace(
        output_create_anchor,
        """\
                % Keep relative paths compatible with the scientific helpers,
                % while placing the actual result beside the PHASE shortcuts.
                outputDir = fullfile('..',outputDir);
                [created,createMessage] = mkdir(outputDir);
                if ~created
                    error('PHASE_Model_beta:outputCreateFailed', ...
                        'Could not create output folder %s: %s',outputDir,createMessage);
                end
                app.outputDir = char(java.io.File(outputDir).getCanonicalPath());
                fprintf('Output folder created: %s\\n',app.outputDir);
""",
        1,
    )

    autoload_start = "            % Automatic load from input_model.mat if it exists\n"
    autoload_end = "        % Value changed function: AOIfiletypeDropDown\n"
    if code.count(autoload_start) != 1 or code.count(autoload_end) != 1:
        raise RuntimeError("Could not locate Model legacy startup-load block")
    start = code.index(autoload_start)
    end = code.index(autoload_end, start)
    code = (
        code[:start]
        + "            % Configuration is applied by the standalone controller.\n"
          "            % Do not push NaN/automatic values through hidden legacy widgets.\n"
          "        end\n\n"
        + code[end:]
    )

    property_anchor = """\
        % Input Files Tab
        filepathIN = ''; % string for .xlsx/.csv path
"""
    if code.count(property_anchor) != 1:
        raise RuntimeError("Could not locate Model backend property anchor")
    code = code.replace(
        property_anchor,
        "        ExternalLogCallback = []\n"
        "        ExternalProgressCallback = []\n"
        "        StopRequested = false\n"
        "        aoi_polygon_lonlat = zeros(0,2)\n\n"
        + property_anchor,
        1,
    )

    methods_anchor = "    % Callbacks that handle component events\n"
    if code.count(methods_anchor) != 1:
        raise RuntimeError("Could not locate Model callback section")
    callback_methods = """\
    methods (Access = public)
        function notifyBetaProgress(app, percentage, phase)
            progress = struct( ...
                'percentage', double(percentage), ...
                'phase', char(string(phase)));
            if ~isempty(app.ExternalProgressCallback)
                try
                    app.ExternalProgressCallback(progress);
                catch
                end
            end
            if ~isempty(app.ExternalLogCallback)
                try
                    app.ExternalLogCallback(progress.phase);
                catch
                end
            end
            drawnow limitrate
            phase_model_beta.throwIfStopped(app);
        end
    end

"""
    code = code.replace(methods_anchor, callback_methods + methods_anchor, 1)

    waitbar_replacements = {
        "            hhh = waitbar(0, 'Preparing environment...');":
            "            app.notifyBetaProgress(0, 'Preparing environment...');",
        "                waitbar(0.1, hhh, 'Loading inputs...');":
            "                app.notifyBetaProgress(10, 'Loading inputs...');",
        "                waitbar(0.2, hhh, 'Importing files...');":
            "                app.notifyBetaProgress(20, 'Importing files...');",
        "                waitbar(0.3, hhh, 'Creating folders...');":
            "                app.notifyBetaProgress(30, 'Creating folders...');",
        "                waitbar(0.4, hhh, 'Preparing data...');":
            "                app.notifyBetaProgress(40, 'Preparing data...');",
        "                waitbar(0.5, hhh, 'Modeling displacement...');":
            "                app.notifyBetaProgress(50, 'Modeling displacement...');",
        "                        waitbar(0.7, hhh, 'Generating report variables...');":
            "                        app.notifyBetaProgress(70, 'Generating report variables...');",
        "                waitbar(0.8, hhh, 'Extrapolating time series...');":
            "                app.notifyBetaProgress(80, 'Extrapolating time series...');",
        "                waitbar(0.9, hhh, 'Generating Excel report...');":
            "                app.notifyBetaProgress(90, 'Generating Excel report...');",
        "                waitbar(1, hhh, 'Processing complete!');":
            "                app.notifyBetaProgress(100, 'Processing complete!');",
    }
    for old, new in waitbar_replacements.items():
        if old not in code:
            raise RuntimeError(f"Could not locate Model progress anchor: {old!r}")
        code = code.replace(old, new)

    code = code.replace(
        "                % close waitbar after a brief delay to show completion\n"
        "                pause(1);\n"
        "                close(hhh);",
        "                % Progress is shown in the PHASE Model run monitor.",
        1,
    )
    error_waitbar_start = (
        "                % 2. check if waitbar exists and is a valid handle before touching it\n"
    )
    error_waitbar_end = "                % 3. update UI label\n"
    if code.count(error_waitbar_start) != 1 or code.count(error_waitbar_end) != 1:
        raise RuntimeError("Could not locate Model error waitbar block")
    start = code.index(error_waitbar_start)
    end = code.index(error_waitbar_end,start)
    code = code[:start] + error_waitbar_end + code[end + len(error_waitbar_end):]
    system_anchor = "[status, cmdout] = system(command);"
    if code.count(system_anchor) != 7:
        raise RuntimeError(
            f"Expected 7 Model report commands, found {code.count(system_anchor)}"
        )
    code = code.replace(
        system_anchor,
        "[status, cmdout] = phase_model_beta.runCommandHidden("
        "app, command, 'Excel report formatting');",
    )
    alert_anchor = (
        "                % 4. show the modal alert with the detailed message\n"
        "                uialert(app.UIFigure, detailedMessage, 'Processing Error', 'Icon', 'error');"
    )
    if code.count(alert_anchor) != 1:
        raise RuntimeError("Could not locate Model processing error handoff")
    code = code.replace(
        alert_anchor,
        "                % 4. hand the error to the visible standalone controller\n"
        "                app.notifyBetaProgress(100, ['Processing failed: ' err.message]);\n"
        "                rethrow(err);",
        1,
    )

    processing_replacements = {
        """\
                    lonMinAOI = app.lonminEditField.Value;
                    lonMaxAOI = app.lonmaxEditField.Value;
                    latMinAOI = app.latminEditField.Value;
                    latMaxAOI = app.latmaxEditField.Value;
                    case false
                    filepathAOI = app.shapefilepathEditField.Value;
""": """\
                    lonMinAOI = app.lonMinAOI;
                    lonMaxAOI = app.lonMaxAOI;
                    latMinAOI = app.latMinAOI;
                    latMaxAOI = app.latMaxAOI;
                    case false
                    filepathAOI = app.filepathAOI;
""",
        """\
                if strcmp(app.MethodDropDown_temporal.Value, 'manual')
                    dtCov_STC1D = app.manualvalueEditField_temporal.Value;
                else
                    dtCov_STC1D = NaN;
                end
""": """\
                dtCov_STC1D = app.dtCov_STC1D;
""",
        """\
                if strcmp(app.MethodDropDown_spatial.Value, 'manual')
                    dsCov_STC1D = app.manualvalueEditField_spatial.Value;
                else
                    dsCov_STC1D = NaN;
                end
""": """\
                dsCov_STC1D = app.dsCov_STC1D;
""",
        "                    varNoise_manual_DET1D = app.manualvalueEditField_noise.Value;":
            "                    varNoise_manual_DET1D = app.varNoise_manual_DET1D;",
        """\
                    num_spl_row_manual_DET1D = app.rownEditField.Value;
                    num_spl_col_manual_DET1D = app.colnEditField.Value;
""": """\
                    num_spl_row_manual_DET1D = app.num_spl_row_manual_DET1D;
                    num_spl_col_manual_DET1D = app.num_spl_col_manual_DET1D;
""",
        "                     lambda_manual_DET1D = app.manualnEditField_4.Value;":
            "                     lambda_manual_DET1D = app.lambda_manual_DET1D;",
        """\
                if strcmp(app.MethodDropDown_temporal_2.Value, 'manual')
                    dtCov_STC2D = app.manualvalueEditField_temporal_2.Value;
                else
                    dtCov_STC2D = NaN;
                end
""": """\
                dtCov_STC2D = app.dtCov_STC2D;
""",
        """\
                if strcmp(app.MethodDropDown_spatial_2.Value, 'manual')
                    dsCov_STC2D = app.manualvalueEditField_spatial_2.Value;
                else
                    dsCov_STC2D = NaN;
                end
""": """\
                dsCov_STC2D = app.dsCov_STC2D;
""",
        "                    varNoise_manual_DET2D = app.manualvalueEditField_noise_2.Value;":
            "                    varNoise_manual_DET2D = app.varNoise_manual_DET2D;",
        """\
                    num_spl_row_manual_DET2D = app.xnEditField.Value;
                    num_spl_col_manual_DET2D = app.ynEditField.Value;
                    num_spl_t_manual_DET2D   = app.tnEditField.Value;
""": """\
                    num_spl_row_manual_DET2D = app.num_spl_row_manual_DET2D;
                    num_spl_col_manual_DET2D = app.num_spl_col_manual_DET2D;
                    num_spl_t_manual_DET2D   = app.num_spl_t_manual_DET2D;
""",
        "                    lambda_manual_DET2D = app.manualnEditField_5.Value;":
            "                    lambda_manual_DET2D = app.lambda_manual_DET2D;",
    }
    for old, new in processing_replacements.items():
        if code.count(old) != 1:
            raise RuntimeError(
                f"Could not locate Model property-processing anchor: {old!r}"
            )
        code = code.replace(old, new, 1)

    polygon_anchor = """\
                    % a) manual coordinates
                    % create the bounding box (polygon) for AOI
                    lonlatAOI = [
                        lonMinAOI, latMinAOI;
                        lonMaxAOI, latMinAOI;
                        lonMaxAOI, latMaxAOI;
                        lonMinAOI, latMaxAOI;
                        lonMinAOI, latMinAOI
                    ];
"""
    polygon_replacement = """\
                    % a) polygon drawn in the standalone map, with the
                    % legacy numeric bounding box retained as fallback
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
"""
    if code.count(polygon_anchor) != 1:
        raise RuntimeError("Could not locate Model AOI polygon anchor")
    code = code.replace(polygon_anchor, polygon_replacement, 1)

    shapefile_start = "                    % b) shapefile\n                    % read bounding box and check for coordinate type\n"
    shapefile_end = "                end\n                \n                % check which PS are inside the AOI\n"
    if code.count(shapefile_start) != 1 or code.count(shapefile_end) != 1:
        raise RuntimeError("Could not locate Model shapefile AOI block")
    start = code.index(shapefile_start)
    end = code.index(shapefile_end,start)
    shapefile_replacement = """\
                    % b) shapefile
                    % Use the same multipart-aware geographic geometry shown
                    % in the standalone map.
                    [lonlatAOI,~,aoiInfo] = phase_model_beta.readAoiShapefile( ...
                        filepathAOI,filepathIN);
                    fprintf('The shapefile AOI coordinates are %s.\\n',aoiInfo.coordinateType);
                    fprintf('AOI polygon parts: %d\\n',aoiInfo.partCount);
                    finiteAOI = all(isfinite(lonlatAOI),2);
                    xyAOI = NaN(size(lonlatAOI));
                    [xAOI,yAOI] = deg2utm( ...
                        lonlatAOI(finiteAOI,2),lonlatAOI(finiteAOI,1));
                    xyAOI(finiteAOI,:) = [xAOI,yAOI];

"""
    code = code[:start] + shapefile_replacement + code[end:]

    shapefile_import_anchor = """\
                % - 2.2) Import the shapefile of the AOI
                if ~flag_AOIbb
                    fileAOI = shaperead(filepathAOI);
                end
"""
    shapefile_import_replacement = """\
                % - 2.2) Import the shapefile of the AOI
                % Loaded later by phase_model_beta.readAoiShapefile so the
                % map and numerical selection share one geometry.
"""
    if code.count(shapefile_import_anchor) != 1:
        raise RuntimeError("Could not locate legacy shapefile import")
    code = code.replace(shapefile_import_anchor,shapefile_import_replacement,1)

    ps_filter_anchor = "                PSidIN_AOI = PSidIN(xyIN_AOI_flag, :);\n"
    if code.count(ps_filter_anchor) != 1:
        raise RuntimeError("Could not locate Model AOI PS-filter anchor")
    code = code.replace(
        ps_filter_anchor,
        ps_filter_anchor
        + "                if isempty(PSidIN_AOI)\n"
          "                    error('PHASE_Model_beta:noPsInsideAoi', ...\n"
          "                        ['The selected AOI contains no persistent scatterers from the ', ...\n"
          "                         'input dataset. Check the AOI shown on the map or select the ', ...\n"
          "                         'full PS extent before starting.']);\n"
          "                end\n",
        1,
    )

    figure_start = "                % - 4.3) Figure of processing scene & AOI\n"
    figure_end = "                fig1_filename = strcat(figsDir, filesep, 'AOI_PS.png');\n"
    if code.count(figure_start) != 1 or code.count(figure_end) != 1:
        raise RuntimeError("Could not locate Model AOI report figure block")
    start = code.index(figure_start)
    end = code.index(figure_end,start)
    figure_replacement = """\
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
"""
    code = code[:start] + figure_replacement + code[end:]

    reverse_anchor = """\
                if ~isnan(ref_centre_lonlat(1))
                    query_lon = ref_centre_lonlat(1);
                    query_lat = ref_centre_lonlat(2);
                else
                    query_lon = mean(lonlatIN_AOI(:,1), 'omitnan');
                    query_lat = mean(lonlatIN_AOI(:,2), 'omitnan');
                end
                \n                % fetch the location data
                [municipality, country] = get_place_from_coordinates(query_lon, query_lat);
"""
    reverse_replacement = """\
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
                    fprintf('Reverse geocoding skipped because the AOI centre is invalid.\\n');
                end
"""
    if code.count(reverse_anchor) != 1:
        raise RuntimeError("Could not locate Model reverse-geocoding block")
    code = code.replace(reverse_anchor,reverse_replacement,1)
    code = code.replace(
        "                aoi_mean_lon = mean(lonlatAOI(:,1));\n"
        "                aoi_mean_lat = mean(lonlatAOI(:,2));",
        "                aoi_mean_lon = mean(lonlatAOI(:,1),'omitnan');\n"
        "                aoi_mean_lat = mean(lonlatAOI(:,2),'omitnan');",
        1,
    )

    grid_start = "                % - 4.4) Create the interpolation grid / centerline based on projDim\n"
    grid_end = "                % - 4.5) Determine municipality and define export filenames\n"
    if code.count(grid_start) != 1 or code.count(grid_end) != 1:
        raise RuntimeError("Could not locate Model interpolation-grid block")
    code = code.replace(
        grid_start,
        grid_start
        + "                if strcmp(procType, 'temporal')\n"
          "                    % Pure temporal modelling works at the observed PS only.\n"
          "                    % Do not allocate an unused centerline or potentially huge 2D grid.\n"
          "                    centerline_data = []; xy_grid = []; lonlat_grid_AOI = [];\n"
          "                    fprintf('Pure temporal mode: spatial grid generation skipped.\\n');\n"
          "                else\n",
        1,
    )
    code = code.replace(grid_end, "                end\n\n" + grid_end, 1)

    figure_exports = {
        "                pause(2)\n                print(f, fig1_filename, '-dpng', '-r300')":
            "                phase_model_beta.exportFigure(f,fig1_filename);",
        "                pause(2)\n                print(f, fig2_filename, '-dpng', '-r300')":
            "                phase_model_beta.exportFigure(f,fig2_filename);",
        "                print(f, fig3_filename, '-dpng', '-r300')":
            "                phase_model_beta.exportFigure(f,fig3_filename);",
        "                    pause(2)\n                    print(f, fig4_filename, '-dpng', '-r300');":
            "                    phase_model_beta.exportFigure(f,fig4_filename);",
        "                    pause(2)\n                    print(f, fig5_filename, '-dpng', '-r300');":
            "                    phase_model_beta.exportFigure(f,fig5_filename);",
    }
    for old, new in figure_exports.items():
        if code.count(old) != 1:
            raise RuntimeError(f"Could not locate report figure export: {old!r}")
        code = code.replace(old,new,1)

    visible_anchor = "            app.UIFigure.Visible = 'on';"
    if code.count(visible_anchor) != 1:
        raise RuntimeError("Could not locate Model engine visibility anchor")
    code = code.replace(
        visible_anchor,
        "            app.UIFigure.Visible = 'off';",
        1,
    )

    threshold_anchor = """\
                    switch coll_proc
                        case 'prediction'
                        OptionalArgs = [OptionalArgs, {'coll_step_est', coll_step_est}];
                    end
"""
    if code.count(threshold_anchor) != 1:
        raise RuntimeError("Could not locate temporal threshold injection anchor")
    threshold_block = threshold_anchor + """\
                    if ismember(procType, {'temporal', 'temporal&NNI'})
                        modelConfig = phase_model_beta.loadConfig( ...
                            phase_model_beta.projectRoot());
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
"""
    code = code.replace(threshold_anchor, threshold_block, 1)
    optional_args_anchor = "                    OptionalArgs = {};"
    if code.count(optional_args_anchor) != 1:
        raise RuntimeError("Could not locate Model temporal OptionalArgs initialisation")
    code = code.replace(
        optional_args_anchor,
        "                    OptionalArgs = {'stop_check', ...\n"
        "                        @() phase_model_beta.throwIfStopped(app)};",
        1,
    )
    # App Designer XML contains pervasive indentation-only/trailing spaces.
    # Normalise them so the editable generated backend remains diff-clean.
    return "\n".join(line.rstrip() for line in code.splitlines()) + "\n"


def main() -> None:
    with zipfile.ZipFile(MLAPP) as archive:
        xml = archive.read("matlab/document.xml").decode("utf-8")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(extract(xml), encoding="utf-8", newline="\n")
    print(f"Generated {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
