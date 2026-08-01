classdef LegacyEngine < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        OutputfolderLabel               matlab.ui.control.Label
        LampLoad                        matlab.ui.control.Lamp
        LampSave                        matlab.ui.control.Lamp
        LoadButton                      matlab.ui.control.Button
        SaveButton                      matlab.ui.control.Button
        StartButton                     matlab.ui.control.Button
        Image_2                         matlab.ui.control.Image
        Image                           matlab.ui.control.Image
        TabGroup                        matlab.ui.container.TabGroup
        InputFilesTab                   matlab.ui.container.Tab
        removealloutputfoldersCheckBox  matlab.ui.control.CheckBox
        pythoninstallationpathEditField  matlab.ui.control.EditField
        pythoninstallationpathEditFieldLabel  matlab.ui.control.Label
        firstdateDatePicker             matlab.ui.control.DatePicker
        firstdateDatePickerLabel        matlab.ui.control.Label
        BrowseButton                    matlab.ui.control.Button
        inputfilepathxlsxcsvEditField   matlab.ui.control.EditField
        inputfilepathxlsxcsvEditFieldLabel  matlab.ui.control.Label
        AOISelectionTab                 matlab.ui.container.Tab
        BrowseShapefileButton           matlab.ui.control.Button
        lonmaxEditField                 matlab.ui.control.NumericEditField
        lonmaxEditFieldLabel            matlab.ui.control.Label
        lonminEditField                 matlab.ui.control.NumericEditField
        lonminEditFieldLabel            matlab.ui.control.Label
        latmaxEditField                 matlab.ui.control.NumericEditField
        latmaxEditFieldLabel            matlab.ui.control.Label
        latminEditField                 matlab.ui.control.NumericEditField
        latminEditFieldLabel            matlab.ui.control.Label
        AOIfiletypeDropDown             matlab.ui.control.DropDown
        AOIfiletypeDropDownLabel        matlab.ui.control.Label
        shapefilepathEditField          matlab.ui.control.EditField
        shapefilepathEditFieldLabel     matlab.ui.control.Label
        ProcessingOptionsTab            matlab.ui.container.Tab
        markersizeforplotsLabel         matlab.ui.control.Label
        EditField                       matlab.ui.control.NumericEditField
        processingtypeDropDown          matlab.ui.control.DropDown
        processingtypeDropDownLabel     matlab.ui.control.Label
        projectdimensionDropDown        matlab.ui.control.DropDown
        projectdimensionDropDownLabel   matlab.ui.control.Label
        AdvancedParametersTab           matlab.ui.container.Tab
        SplinesPanel_spatial_2D         matlab.ui.container.Panel
        tnEditField                     matlab.ui.control.NumericEditField
        tnEditFieldLabel                matlab.ui.control.Label
        ynEditField                     matlab.ui.control.NumericEditField
        ynEditFieldLabel                matlab.ui.control.Label
        manualvalueEditField_noise_2    matlab.ui.control.NumericEditField
        manualvalueEditField_3Label_4   matlab.ui.control.Label
        NoiseDropDown_2                 matlab.ui.control.DropDown
        NoiseDropDown_2Label            matlab.ui.control.Label
        manualnEditField_5              matlab.ui.control.NumericEditField
        manualnEditField_5Label         matlab.ui.control.Label
        LambdaDropDown_3                matlab.ui.control.DropDown
        LambdaDropDown_3Label           matlab.ui.control.Label
        IndexDropDown_3                 matlab.ui.control.DropDown
        IndexDropDown_3Label            matlab.ui.control.Label
        xnEditField                     matlab.ui.control.NumericEditField
        xnEditFieldLabel                matlab.ui.control.Label
        MethodDropDown_5                matlab.ui.control.DropDown
        MethodDropDown_5Label_5         matlab.ui.control.Label
        CovariancePanel_2D              matlab.ui.container.Panel
        ModelDropDown_spatial_2         matlab.ui.control.DropDown
        ModelDropDownLabel_4            matlab.ui.control.Label
        ModelDropDown_temporal_2        matlab.ui.control.DropDown
        ModelDropDownLabel_3            matlab.ui.control.Label
        manualvalueEditField_spatial_2  matlab.ui.control.NumericEditField
        manualvalueEditField_3Label_3   matlab.ui.control.Label
        MethodDropDown_spatial_2        matlab.ui.control.DropDown
        MethodDropDown_5Label_3         matlab.ui.control.Label
        SpatialLabel_2                  matlab.ui.control.Label
        TemporalLabel_2                 matlab.ui.control.Label
        manualvalueEditField_temporal_2  matlab.ui.control.NumericEditField
        manualvalueEditField_2Label_2   matlab.ui.control.Label
        MethodDropDown_temporal_2       matlab.ui.control.DropDown
        MethodDropDown_4Label_3         matlab.ui.control.Label
        SplinesPanel_spatial            matlab.ui.container.Panel
        colnEditField                   matlab.ui.control.NumericEditField
        colnEditFieldLabel              matlab.ui.control.Label
        manualvalueEditField_noise      matlab.ui.control.NumericEditField
        manualvalueEditField_3Label_2   matlab.ui.control.Label
        NoiseDropDown                   matlab.ui.control.DropDown
        MethodDropDown_5Label_2         matlab.ui.control.Label
        manualnEditField_4              matlab.ui.control.NumericEditField
        manualnEditField_4Label         matlab.ui.control.Label
        LambdaDropDown_2                matlab.ui.control.DropDown
        LambdaDropDown_2Label           matlab.ui.control.Label
        IndexDropDown_2                 matlab.ui.control.DropDown
        IndexDropDown_2Label            matlab.ui.control.Label
        rownEditField                   matlab.ui.control.NumericEditField
        rownEditFieldLabel              matlab.ui.control.Label
        MethodDropDown_4                matlab.ui.control.DropDown
        MethodDropDown_4Label_2         matlab.ui.control.Label
        CovariancePanel                 matlab.ui.container.Panel
        ModelDropDown_spatial           matlab.ui.control.DropDown
        ModelDropDownLabel_2            matlab.ui.control.Label
        ModelDropDown_temporal          matlab.ui.control.DropDown
        ModelDropDownLabel              matlab.ui.control.Label
        manualvalueEditField_spatial    matlab.ui.control.NumericEditField
        manualvalueEditField_3Label     matlab.ui.control.Label
        MethodDropDown_spatial          matlab.ui.control.DropDown
        MethodDropDown_5Label           matlab.ui.control.Label
        SpatialLabel                    matlab.ui.control.Label
        TemporalLabel                   matlab.ui.control.Label
        manualvalueEditField_temporal   matlab.ui.control.NumericEditField
        manualvalueEditField_2Label     matlab.ui.control.Label
        MethodDropDown_temporal         matlab.ui.control.DropDown
        MethodDropDown_4Label           matlab.ui.control.Label
        GeneralParametersPanel          matlab.ui.container.Panel
        TimestepEditField               matlab.ui.control.NumericEditField
        TimestepEditFieldLabel          matlab.ui.control.Label
        TiltedepochwiseplanesCheckBox   matlab.ui.control.CheckBox
        PolynomialmaxdegreeDropDown     matlab.ui.control.DropDown
        PolynomialmaxdegreeDropDownLabel  matlab.ui.control.Label
        ResidualatmosphericartifactsDropDown  matlab.ui.control.DropDown
        ResidualatmosphericartifactsDropDownLabel  matlab.ui.control.Label
        GridresolutionmEditField        matlab.ui.control.NumericEditField
        GridresolutionmEditFieldLabel   matlab.ui.control.Label
        CenterlineresolutionmEditField  matlab.ui.control.NumericEditField
        CenterlineresolutionmEditFieldLabel  matlab.ui.control.Label
        MinperiodmonthsforoutlierdetectionwithsplinesEditField  matlab.ui.control.NumericEditField
        MinperiodmonthsforoutlierdetectionwithsplinesEditFieldLabel  matlab.ui.control.Label
        CollocationPanel                matlab.ui.container.Panel
        estimationstepEditField         matlab.ui.control.NumericEditField
        estimationstepEditFieldLabel    matlab.ui.control.Label
        MethodDropDown_3                matlab.ui.control.DropDown
        MethodDropDown_3Label           matlab.ui.control.Label
        SplinesPanel                    matlab.ui.container.Panel
        manualnEditField_2              matlab.ui.control.NumericEditField
        manualnEditField_2Label         matlab.ui.control.Label
        LambdaDropDown                  matlab.ui.control.DropDown
        LambdaDropDownLabel             matlab.ui.control.Label
        IndexDropDown                   matlab.ui.control.DropDown
        IndexDropDownLabel              matlab.ui.control.Label
        manualnEditField                matlab.ui.control.NumericEditField
        manualnEditFieldLabel           matlab.ui.control.Label
        MethodDropDown_2                matlab.ui.control.DropDown
        MethodDropDown_2Label           matlab.ui.control.Label
        NoiseVariancePanel              matlab.ui.container.Panel
        nlooksEditField                 matlab.ui.control.NumericEditField
        nlooksEditFieldLabel            matlab.ui.control.Label
        constellationDropDown           matlab.ui.control.DropDown
        constellationDropDownLabel      matlab.ui.control.Label
        BrowseButton_2                  matlab.ui.control.Button
        cohfolderEditField              matlab.ui.control.EditField
        cohfolderEditFieldLabel         matlab.ui.control.Label
        manualvalueEditField            matlab.ui.control.NumericEditField
        manualvalueEditFieldLabel       matlab.ui.control.Label
        MethodDropDown                  matlab.ui.control.DropDown
        MethodDropDownLabel             matlab.ui.control.Label
        QuerypointsTab                  matlab.ui.container.Tab
        exportobsCheckBox               matlab.ui.control.CheckBox
        BrowseExtrapolationButton       matlab.ui.control.Button
        timeseriesinterpolationatquerypointsCheckBox  matlab.ui.control.CheckBox
        interpolationfilepathEditField  matlab.ui.control.EditField
        interpolationfilepathEditFieldLabel  matlab.ui.control.Label
    end


    properties (Access = public)

        ExternalLogCallback = []
        ExternalProgressCallback = []
        StopRequested = false
        aoi_polygon_lonlat = zeros(0,2)

        % Input Files Tab
        filepathIN = ''; % string for .xlsx/.csv path
        pythonPath = ''; % string for Python path
        t0IN = datetime('now'); % datetime for first date
        flag_ouputDir = false; % logical for removing output folders
        outputDir % output folder name (e.g., 'output_001')

        % AOI Selection Tab
        flag_AOIbb = false; % logical for Shapefile vs Bounding Box
        filepathAOI = ''; % string for shapefile path
        lonMinAOI = []; % numeric for bounding box
        lonMaxAOI = []; % numeric
        latMinAOI = []; % numeric
        latMaxAOI = []; % numeric

        % Processing Options Tab
        projDim = '1D'; % string: '1D' or '2D'
        procType = 'temporal'; % string: 'temporal', 'temporal&NNI', 'spatialDET', 'spatialSTC'

        % Extrapolation Tab
        flag_PSinterp = false; % logical for PS interpolation
        flag_tsExtr = false; % logical for extrapolation
        filepath_EXTR = ''; % string for extrapolation file path

        % Advanced Parameters
        varNoise_method = 'auto'; % string: 'auto', 'manual', 'coherence'
        varNoise_manual = NaN; % numeric, NaN if not manual
        coherence_dir = ''; % string, empty if not coherence
        constellation = 'Sentinel-1'; % string
        num_looks = []; % numeric, empty if not coherence
        num_spl_method = 'auto'; % string: 'auto', 'manual'
        spline_method = 'variance'; % string: 'variance', 'MDL', 'F_test', 'chi2_test'
        num_spl_manual = NaN; % numeric, NaN if not manual
        lambda_method = 'auto'; % string: 'auto', 'manual'
        lambda_manual = NaN; % numeric, NaN if not manual
        coll_proc = 'filtering'; % string: 'filtering', 'prediction'
        coll_step_est = NaN; % numeric, NaN if not prediction
        minMonths = []; % numeric
        cline_resolution = []; % numeric, for 1D
        grid_resolution = []; % numeric, for 2D
        dtCov_STC1D = NaN; % numeric, NaN for auto
        dsCov_STC1D = NaN; % numeric, NaN for auto
        varNoise_DET1D = 'auto'; % string: 'auto', 'manual'
        varNoise_manual_DET1D = NaN; % numeric, NaN if not manual
        num_spl_method_DET1D = 'auto'; % string: 'auto', 'manual'
        spline_method_DET1D = 'MDL'; % string: 'variance', 'MDL', 'F_test', 'chi2_test'
        num_spl_row_manual_DET1D = NaN; % numeric, NaN if not manual
        num_spl_col_manual_DET1D = NaN; % numeric, NaN if not manual
        lambda_method_DET1D = 'auto'; % string: 'auto', 'manual'
        lambda_manual_DET1D = NaN; % numeric, NaN if not manual
        step_t_ST = 1; % numeric
        dtCov_STC2D = NaN; % numeric, NaN for auto
        dsCov_STC2D = NaN; % numeric, NaN for auto
        varNoise_DET2D = 'auto'; % string: 'auto', 'manual'
        varNoise_manual_DET2D = NaN; % numeric, NaN if not manual
        num_spl_method_DET2D = 'auto'; % string: 'auto', 'manual'
        spline_method_DET2D = 'MDL'; % string: 'variance', 'MDL', 'F_test', 'chi2_test'
        num_spl_row_manual_DET2D = NaN; % numeric, NaN if not manual
        num_spl_col_manual_DET2D = NaN; % numeric, NaN if not manual
        num_spl_t_manual_DET2D = NaN; % numeric, NaN if not manual
        lambda_method_DET2D = 'auto'; % string: 'auto', 'manual'
        lambda_manual_DET2D = NaN; % numeric, NaN if not manual
        detrendMethod = 'no'; % string: 'yes', 'no'
        markerSize = 10; % Default value for marker size

        % Covariance Models defaults
        tCovModel_STC1D = 'gaussian';
        sCovModel_STC1D = 'gaussian';
        tCovModel_STC2D = 'gaussian';
        sCovModel_STC2D = 'gaussian';

        % Advanced Options
        polyDegreeST = '1'; % Default string for DropDown
        useInclinedMeansST = true; % Logical for checkbox

    end



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

    % Callbacks that handle component events
    methods (Access = public)

        % Code that executes after component creation
        function startupFcn(app)

            rootDir = phase_model_beta.projectRoot();
            cd(rootDir);
            addpath(fullfile(rootDir, 'MatlabFunctions'));
            phase_model_beta.themeLegacyEngine(app);
            app.pythonPath = phase_model_beta.resolvePythonPath(app.pythonPath);
            app.pythoninstallationpathEditField.Value = app.pythonPath;

            % Initialize visibility
            app.shapefilepathEditField.Visible = 'on';
            app.shapefilepathEditFieldLabel.Visible = 'on';
            app.BrowseShapefileButton.Visible = 'on';
            app.lonminEditField.Visible = 'off';
            app.lonminEditFieldLabel.Visible = 'off';
            app.lonmaxEditField.Visible = 'off';
            app.lonmaxEditFieldLabel.Visible = 'off';
            app.latminEditField.Visible = 'off';
            app.latminEditFieldLabel.Visible = 'off';
            app.latmaxEditField.Visible = 'off';
            app.latmaxEditFieldLabel.Visible = 'off';
            app.interpolationfilepathEditField.Visible = 'off';
            app.interpolationfilepathEditFieldLabel.Visible = 'off';
            app.BrowseExtrapolationButton.Visible = 'off';
            % Advanced Parameters visibility
            app.manualvalueEditField.Visible = 'off';
            app.manualvalueEditFieldLabel.Visible = 'off';
            app.cohfolderEditField.Visible = 'off';
            app.cohfolderEditFieldLabel.Visible = 'off';
            app.BrowseButton_2.Visible = 'off';
            app.constellationDropDown.Visible = 'off';
            app.constellationDropDownLabel.Visible = 'off';
            app.nlooksEditField.Visible = 'off';
            app.nlooksEditFieldLabel.Visible = 'off';
            app.manualnEditField.Visible = 'off';
            app.manualnEditFieldLabel.Visible = 'off';
            app.IndexDropDown.Visible = 'on';
            app.IndexDropDownLabel.Visible = 'on';
            app.manualnEditField_2.Visible = 'off';
            app.manualnEditField_2Label.Visible = 'off';
            app.estimationstepEditField.Visible = 'off';
            app.estimationstepEditFieldLabel.Visible = 'off';
            app.GridresolutionmEditField.Visible = 'off';
            app.GridresolutionmEditFieldLabel.Visible = 'off';
            app.manualvalueEditField_temporal.Visible = 'off';
            app.manualvalueEditField_2Label.Visible = 'off';
            app.manualvalueEditField_temporal_2.Visible = 'off';
            app.manualvalueEditField_2Label_2.Visible = 'off';
            app.manualvalueEditField_spatial.Visible = 'off';
            app.manualvalueEditField_3Label.Visible = 'off';
            app.manualvalueEditField_spatial_2.Visible = 'off';
            app.manualvalueEditField_3Label_3.Visible = 'off';
            app.rownEditField.Visible = 'off';
            app.rownEditFieldLabel.Visible = 'off';
            app.xnEditField.Visible = 'off';
            app.xnEditFieldLabel.Visible = 'off';
            app.colnEditField.Visible = 'off';
            app.colnEditFieldLabel.Visible = 'off';
            app.ynEditField.Visible = 'off';
            app.ynEditFieldLabel.Visible = 'off';
            app.tnEditField.Visible = 'off';
            app.tnEditFieldLabel.Visible = 'off';
            app.manualnEditField_4.Visible = 'off';
            app.manualnEditField_4Label.Visible = 'off';
            app.manualnEditField_5.Visible = 'off';
            app.manualnEditField_5Label.Visible = 'off';
            app.manualvalueEditField_noise.Visible = 'off';
            app.manualvalueEditField_noise_2.Visible = 'off';
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Visible = 'off';
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditFieldLabel.Visible = 'off';
            app.manualvalueEditField_3Label_2.Visible = 'off';
            app.IndexDropDown_2.Visible = 'off';
            app.IndexDropDown_2Label.Visible = 'off';
            app.IndexDropDown_3.Visible = 'off';
            app.IndexDropDown_3Label.Visible = 'off';
            app.TimestepEditField.Visible = 'off';
            app.TimestepEditFieldLabel.Visible = 'off';
            app.manualvalueEditField_3Label_4.Visible = 'off';
            app.ResidualatmosphericartifactsDropDown.Visible = 'off';
            app.ResidualatmosphericartifactsDropDownLabel.Visible = 'off';
            app.TiltedepochwiseplanesCheckBox.Visible = 'off';
            app.PolynomialmaxdegreeDropDown.Visible = 'off';
            app.PolynomialmaxdegreeDropDownLabel.Visible = 'off';
            % Panel visibility
            app.NoiseVariancePanel.Visible = 'on';
            app.SplinesPanel.Visible = 'on';
            app.CollocationPanel.Visible = 'on';
            app.GeneralParametersPanel.Visible = 'on';
            app.CovariancePanel.Visible = 'off';
            app.SplinesPanel_spatial.Visible = 'off';
            app.CovariancePanel_2D.Visible = 'off';
            app.SplinesPanel_spatial_2D.Visible = 'off';

            % Configuration is applied by the standalone controller.
            % Do not push NaN/automatic values through hidden legacy widgets.
        end

        % Value changed function: AOIfiletypeDropDown
        function AOIfiletypeDropDownValueChanged(app, event)
            if strcmp(app.AOIfiletypeDropDown.Value, 'Shapefile')
                % Show shapefile fields
                app.shapefilepathEditField.Visible = 'on';
                app.shapefilepathEditFieldLabel.Visible = 'on';
                app.BrowseShapefileButton.Visible = 'on';
                % Hide bounding box fields
                app.lonminEditField.Visible = 'off';
                app.lonminEditFieldLabel.Visible = 'off';
                app.lonmaxEditField.Visible = 'off';
                app.lonmaxEditFieldLabel.Visible = 'off';
                app.latminEditField.Visible = 'off';
                app.latminEditFieldLabel.Visible = 'off';
                app.latmaxEditField.Visible = 'off';
                app.latmaxEditFieldLabel.Visible = 'off';
                app.flag_AOIbb = false;
            else % Bounding Box
                % Hide shapefile fields
                app.shapefilepathEditField.Visible = 'off';
                app.shapefilepathEditFieldLabel.Visible = 'off';
                app.BrowseShapefileButton.Visible = 'off';
                % Show bounding box fields
                app.lonminEditField.Visible = 'on';
                app.lonminEditFieldLabel.Visible = 'on';
                app.lonmaxEditField.Visible = 'on';
                app.lonmaxEditFieldLabel.Visible = 'on';
                app.latminEditField.Visible = 'on';
                app.latminEditFieldLabel.Visible = 'on';
                app.latmaxEditField.Visible = 'on';
                app.latmaxEditFieldLabel.Visible = 'on';
                app.flag_AOIbb = true;
            end
        end

        % Button pushed function: BrowseShapefileButton
        function BrowseShapefileButtonPushed(app, event)
            [file, path] = uigetfile('*.shp', 'Select Shapefile');
            if ~isequal(file, 0)
                app.shapefilepathEditField.Value = fullfile(path, file);
                app.filepathAOI = app.shapefilepathEditField.Value;
            end
        end

        % Value changed function:
        % timeseriesinterpolationatquerypointsCheckBox
        function timeseriesinterpolationatquerypointsCheckBoxValueChanged(app, event)
            app.flag_tsExtr = app.timeseriesinterpolationatquerypointsCheckBox.Value;
            if app.flag_tsExtr
                % Show extrapolation fields
                app.interpolationfilepathEditField.Visible = 'on';
                app.interpolationfilepathEditFieldLabel.Visible = 'on';
                app.BrowseExtrapolationButton.Visible = 'on';
            else
                % Hide extrapolation fields
                app.interpolationfilepathEditField.Visible = 'off';
                app.interpolationfilepathEditFieldLabel.Visible = 'off';
                app.BrowseExtrapolationButton.Visible = 'off';
            end
        end

        % Button pushed function: BrowseExtrapolationButton
        function BrowseExtrapolationButtonPushed(app, event)
            [file, path] = uigetfile('*.txt', 'Select Interpolation File');
            if ~isequal(file, 0)
                app.interpolationfilepathEditField.Value = fullfile(path, file);
                app.filepath_EXTR = app.interpolationfilepathEditField.Value;
            end
        end

        % Button pushed function: BrowseButton
        function BrowseButtonPushed(app, event)
            [file, path] = uigetfile({'*.xlsx;*.csv'}, 'Select Time Series File');
            if ~isequal(file, 0)
                app.inputfilepathxlsxcsvEditField.Value = fullfile(path, file);
                app.filepathIN = app.inputfilepathxlsxcsvEditField.Value;
            end
        end

        % Value changed function: inputfilepathxlsxcsvEditField
        function inputfilepathxlsxcsvEditFieldValueChanged(app, event)
            app.filepathIN = app.inputfilepathxlsxcsvEditField.Value;
        end

        % Value changed function: firstdateDatePicker
        function firstdateDatePickerValueChanged(app, event)
            app.t0IN = app.firstdateDatePicker.Value;
        end

        % Value changed function: pythoninstallationpathEditField
        function pythoninstallationpathEditFieldValueChanged(app, event)
            app.pythonPath = app.pythoninstallationpathEditField.Value;
        end

        % Value changed function: removealloutputfoldersCheckBox
        function removealloutputfoldersCheckBoxValueChanged(app, event)
            app.flag_ouputDir = app.removealloutputfoldersCheckBox.Value;
        end

        % Value changed function: projectdimensionDropDown
        function projectdimensionDropDownValueChanged(app, event)
            app.projDim = app.projectdimensionDropDown.Value;
            app.GridresolutionmEditField.Visible = strcmp(app.projDim, '2D');
            app.GridresolutionmEditFieldLabel.Visible = strcmp(app.projDim, '2D');
            app.CenterlineresolutionmEditField.Visible = strcmp(app.projDim, '1D');
            app.CenterlineresolutionmEditFieldLabel.Visible = strcmp(app.projDim, '1D');
            if strcmp(app.projDim, '1D')
                app.grid_resolution = NaN;
                app.GridresolutionmEditField.Value = 0; % Clear UI field
            else
                app.cline_resolution = NaN;
                app.CenterlineresolutionmEditField.Value = 0; % Clear UI field
            end
            app.processingtypeDropDownValueChanged([]); % Update panel layout
        end

        % Value changed function: processingtypeDropDown
        function processingtypeDropDownValueChanged(app, event)

            switch app.processingtypeDropDown.Value
                case 'Temporal'
                    app.procType = 'temporal';
                case 'Temporal & NNI'
                    app.procType = 'temporal&NNI';
                case 'Spatio-temporal Deterministic'
                    app.procType = 'spatialDET';
                case 'Spatio-temporal Stochastic'
                    app.procType = 'spatialSTC';
            end

            % Determine panels to show based on procType
            switch app.procType
                case {'temporal', 'temporal&NNI'}
                    panelsToShow = {'NoiseVariancePanel', 'SplinesPanel', 'CollocationPanel', 'GeneralParametersPanel'};
                case 'spatialDET'
                    if strcmp(app.projDim, '1D')
                        panelsToShow = {'SplinesPanel_spatial', 'GeneralParametersPanel'};
                    elseif strcmp(app.projDim, '2D')
                        panelsToShow = {'SplinesPanel_spatial_2D', 'GeneralParametersPanel'};
                    end
                case 'spatialSTC'
                    if strcmp(app.projDim, '1D')
                        panelsToShow = {'CovariancePanel', 'GeneralParametersPanel'};
                    elseif strcmp(app.projDim, '2D')
                        panelsToShow = {'CovariancePanel_2D', 'GeneralParametersPanel'};
                    end
            end

            % Hide all panels
            app.NoiseVariancePanel.Visible = 'off';
            app.SplinesPanel.Visible = 'off';
            app.CollocationPanel.Visible = 'off';
            app.GeneralParametersPanel.Visible = 'off';
            app.CovariancePanel.Visible = 'off';
            app.SplinesPanel_spatial.Visible = 'off';
            app.CovariancePanel_2D.Visible = 'off';
            app.SplinesPanel_spatial_2D.Visible = 'off';

            % Assign positions based on number of panels
            margin = 30;
            panelWidth = 460;
            panelHeight = 180;
            positions = [
                margin, 454 - margin - panelHeight, panelWidth, panelHeight; % (1,1) top-left
                1000 - margin - panelWidth, 454 - margin - panelHeight, panelWidth, panelHeight; % (1,2) top-right
                margin, margin, panelWidth, panelHeight; % (2,1) bottom-left
                1000 - margin - panelWidth, margin, panelWidth, panelHeight % (2,2) bottom-right
            ];
            numPanels = length(panelsToShow);
            for i = 1:numPanels
                app.(panelsToShow{i}).Visible = 'on';
                app.(panelsToShow{i}).Position = positions(i, :);
            end

            % Set NaN for unused properties and clear UI fields
            if ~strcmp(app.procType, 'spatialSTC') && strcmp(app.projDim, '1D')
                app.dtCov_STC1D = NaN;
                app.dsCov_STC1D = NaN;
                app.manualvalueEditField_temporal.Value = 0; % Clear UI field
                app.manualvalueEditField_spatial.Value = 0; % Clear UI field
                app.tCovModel_STC1D = 'gaussian';
                app.sCovModel_STC1D = 'gaussian';
            end
            if ~strcmp(app.procType, 'spatialDET') && strcmp(app.projDim, '1D')
                app.varNoise_manual_DET1D = NaN;
                app.num_spl_row_manual_DET1D = NaN;
                app.num_spl_col_manual_DET1D = NaN;
                app.lambda_manual_DET1D = NaN;
                app.manualvalueEditField_noise.Value = 0; % Clear UI field
                app.rownEditField.Value = 2; % Default to minimum valid value
                app.colnEditField.Value = 2; % Default to minimum valid value
                app.manualnEditField_4.Value = 0; % Clear UI field
            end
            if ~strcmp(app.procType, 'spatialSTC') && strcmp(app.projDim, '2D')
                app.dtCov_STC2D = NaN;
                app.dsCov_STC2D = NaN;
                app.manualvalueEditField_temporal_2.Value = 0; % Clear UI field
                app.manualvalueEditField_spatial_2.Value = 0; % Clear UI field
                app.tCovModel_STC2D = 'gaussian';
                app.sCovModel_STC2D = 'gaussian';
            end
            if ~strcmp(app.procType, 'spatialDET') && strcmp(app.projDim, '2D')
                app.varNoise_manual_DET2D = NaN;
                app.num_spl_row_manual_DET2D = NaN;
                app.num_spl_col_manual_DET2D = NaN;
                app.num_spl_t_manual_DET2D = NaN;
                app.lambda_manual_DET2D = NaN;
                app.manualvalueEditField_noise_2.Value = 0; % Clear UI field
                app.xnEditField.Value = 2; % Default to minimum valid value
                app.ynEditField.Value = 2; % Default to minimum valid value
                app.tnEditField.Value = 2; % Default to minimum valid value
                app.manualnEditField_5.Value = 0; % Clear UI field
            end
            % Control visibility of MinperiodmonthsforoutlierdetectionwithsplinesEditField
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditFieldLabel.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
            if ismember(app.procType, {'temporal', 'temporal&NNI'})
                app.minMonths = NaN; % Clear internal property
                app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Value = 0; % Clear UI field
            end
            % Control visibility for additional general ST parameters
            app.TimestepEditField.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
            app.TimestepEditFieldLabel.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
            if ismember(app.procType, {'temporal', 'temporal&NNI'})
                app.step_t_ST = 1; % Clear internal property
                app.TimestepEditField.Value = 1; % Clear UI field
            end

            % Handle Atmospheric Artifacts Dropdown Visibility
            % Hide for Temporal, Show for Spatial
            isTemporal = ismember(app.procType, {'temporal', 'temporal&NNI'});
            app.ResidualatmosphericartifactsDropDown.Visible = ~isTemporal;
            app.ResidualatmosphericartifactsDropDownLabel.Visible = ~isTemporal;

            if isTemporal
                % Force hide everything if in Temporal mode
                app.TiltedepochwiseplanesCheckBox.Visible = 'off';
                app.PolynomialmaxdegreeDropDown.Visible = 'off';
                app.PolynomialmaxdegreeDropDownLabel.Visible = 'off';

                % Reset to default "cleanObs" (No artifact removal) internally
                app.detrendMethod = 'cleanObs';
                app.ResidualatmosphericartifactsDropDown.Value = 'no';
            else
                % Only reset to defaults if a human manually triggered the change
                if ~isempty(event)
                    app.ResidualatmosphericartifactsDropDown.Value = 'no';
                    app.detrendMethod = 'residualAtm';
                    app.useInclinedMeansST = true;
                    app.TiltedepochwiseplanesCheckBox.Value = false;
                    app.polyDegreeST = '1';
                    app.PolynomialmaxdegreeDropDown.Value = '1';
                    % Spatial Mode -> Trigger the Dropdown Logic to decide what to show
                    % Call the other callback manually to enforce the logic based on the current dropdown value
                    app.ResidualatmosphericartifactsDropDownValueChanged([]);
                        % Set internal defaults if not already set
                        if isempty(app.tCovModel_STC1D), app.tCovModel_STC1D = 'gaussian'; end
                        if isempty(app.sCovModel_STC1D), app.sCovModel_STC1D = 'gaussian'; end
                        if isempty(app.tCovModel_STC2D), app.tCovModel_STC2D = 'gaussian'; end
                        if isempty(app.sCovModel_STC2D), app.sCovModel_STC2D = 'gaussian'; end

                        % Safe assignment for Temporal Model 1D
                        if ismember(app.tCovModel_STC1D, app.ModelDropDown_temporal.Items)
                            app.ModelDropDown_temporal.Value = app.tCovModel_STC1D;
                        else
                            % Default to the first item if the stored value is invalid
                            app.ModelDropDown_temporal.Value = app.ModelDropDown_temporal.Items{1};
                            app.tCovModel_STC1D = app.ModelDropDown_temporal.Items{1}; % Update internal var
                        end

                        % Safe assignment for Spatial Model 1D
                        if ismember(app.sCovModel_STC1D, app.ModelDropDown_spatial.Items)
                            app.ModelDropDown_spatial.Value = app.sCovModel_STC1D;
                        else
                            app.ModelDropDown_spatial.Value = app.ModelDropDown_spatial.Items{1};
                            app.sCovModel_STC1D = app.ModelDropDown_spatial.Items{1};
                        end

                        % Safe assignment for Temporal Model 2D
                        if ismember(app.tCovModel_STC2D, app.ModelDropDown_temporal_2.Items)
                            app.ModelDropDown_temporal_2.Value = app.tCovModel_STC2D;
                        else
                            app.ModelDropDown_temporal_2.Value = app.ModelDropDown_temporal_2.Items{1};
                            app.tCovModel_STC2D = app.ModelDropDown_temporal_2.Items{1};
                        end

                        % Safe assignment for Spatial Model 2D
                        if ismember(app.sCovModel_STC2D, app.ModelDropDown_spatial_2.Items)
                            app.ModelDropDown_spatial_2.Value = app.sCovModel_STC2D;
                        else
                            app.ModelDropDown_spatial_2.Value = app.ModelDropDown_spatial_2.Items{1};
                            app.sCovModel_STC2D = app.ModelDropDown_spatial_2.Items{1};
                        end
                    app.ModelDropDown_temporal.Value = app.tCovModel_STC1D;
                    app.ModelDropDown_spatial.Value = app.sCovModel_STC1D;
                    app.ModelDropDown_temporal_2.Value = app.tCovModel_STC2D;
                    app.ModelDropDown_spatial_2.Value = app.sCovModel_STC2D;
                end
                % Set internal defaults if not already set
                if isempty(app.tCovModel_STC1D), app.tCovModel_STC1D = 'gaussian'; end
            end

        end

        % Value changed function: EditField
        function EditFieldValueChanged(app, event)
            app.markerSize = app.EditField.Value;
        end

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
            config = struct();
            config.filepathIN = app.filepathIN;
            config.pythonPath = app.pythonPath;
            config.t0IN = app.t0IN;
            config.flag_ouputDir = app.flag_ouputDir;
            config.flag_AOIbb = app.flag_AOIbb;
            config.filepathAOI = app.filepathAOI;
            config.lonMinAOI = app.lonMinAOI;
            config.lonMaxAOI = app.lonMaxAOI;
            config.latMinAOI = app.latMinAOI;
            config.latMaxAOI = app.latMaxAOI;
            config.projDim = app.projDim;
            config.procType = app.procType;
            config.markerSize = app.markerSize;
            config.flag_PSinterp = app.flag_PSinterp;
            config.flag_tsExtr = app.flag_tsExtr;
            config.filepath_EXTR = app.filepath_EXTR;
            config.varNoise_method = app.varNoise_method;
            config.varNoise_manual = app.varNoise_manual;
            config.coherence_dir = app.coherence_dir;
            config.constellation = app.constellation;
            config.num_looks = app.num_looks;
            config.num_spl_method = app.num_spl_method;
            config.spline_method = app.spline_method;
            config.num_spl_manual = app.num_spl_manual;
            config.lambda_method = app.lambda_method;
            config.lambda_manual = app.lambda_manual;
            config.coll_proc = app.coll_proc;
            config.coll_step_est = app.coll_step_est;
            config.minMonths = app.minMonths;
            config.cline_resolution = app.cline_resolution;
            config.grid_resolution = app.grid_resolution;
            config.dtCov_STC1D = app.dtCov_STC1D;
            config.dsCov_STC1D = app.dsCov_STC1D;
            config.tCovModel_STC2D = app.tCovModel_STC2D;
            config.tCovModel_STC1D = app.tCovModel_STC1D;
            config.sCovModel_STC2D = app.sCovModel_STC2D;
            config.sCovModel_STC1D = app.sCovModel_STC1D;
            config.varNoise_DET1D = app.varNoise_DET1D;
            config.varNoise_manual_DET1D = app.varNoise_manual_DET1D;
            config.num_spl_method_DET1D = app.num_spl_method_DET1D;
            config.spline_method_DET1D = app.spline_method_DET1D;
            config.num_spl_row_manual_DET1D = app.num_spl_row_manual_DET1D;
            config.num_spl_col_manual_DET1D = app.num_spl_col_manual_DET1D;
            config.lambda_method_DET1D = app.lambda_method_DET1D;
            config.lambda_manual_DET1D = app.lambda_manual_DET1D;
            config.step_t_ST = app.step_t_ST;
            config.dtCov_STC2D = app.dtCov_STC2D;
            config.dsCov_STC2D = app.dsCov_STC2D;
            config.varNoise_DET2D = app.varNoise_DET2D;
            config.varNoise_manual_DET2D = app.varNoise_manual_DET2D;
            config.num_spl_method_DET2D = app.num_spl_method_DET2D;
            config.spline_method_DET2D = app.spline_method_DET2D;
            config.num_spl_row_manual_DET2D = app.num_spl_row_manual_DET2D;
            config.num_spl_col_manual_DET2D = app.num_spl_col_manual_DET2D;
            config.num_spl_t_manual_DET2D = app.num_spl_t_manual_DET2D;
            config.lambda_method_DET2D = app.lambda_method_DET2D;
            config.lambda_manual_DET2D = app.lambda_manual_DET2D;
            config.useInclinedMeansST = app.useInclinedMeansST;
            config.markerSize = app.markerSize;
            config.polyDegreeST = app.polyDegreeST;
            config.detrendMethod = app.detrendMethod;
            appDir = phase_model_beta.projectRoot();
            matFile = fullfile(appDir, 'input_model.mat');
            try
                save(matFile, 'config');
                app.LampSave.Color = 'green';
            catch
                app.LampSave.Color = 'red';
            end
        end

        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
            appDir = phase_model_beta.projectRoot();
            matFile = fullfile(appDir, 'input_model.mat');
            if exist(matFile, 'file') == 2
                try
                    load(matFile, 'config');
                    app.filepathIN = config.filepathIN;
                    app.pythonPath = config.pythonPath;
                    app.t0IN = config.t0IN;
                    app.flag_ouputDir = config.flag_ouputDir;
                    app.flag_AOIbb = config.flag_AOIbb;
                    app.filepathAOI = config.filepathAOI;
                    app.lonMinAOI = config.lonMinAOI;
                    app.lonMaxAOI = config.lonMaxAOI;
                    app.latMinAOI = config.latMinAOI;
                    app.latMaxAOI = config.latMaxAOI;
                    app.projDim = config.projDim;
                    app.procType = config.procType;
                    app.markerSize = config.markerSize;
                    app.flag_PSinterp = config.flag_PSinterp;
                    app.flag_tsExtr = config.flag_tsExtr;
                    app.filepath_EXTR = config.filepath_EXTR;
                    app.varNoise_method = config.varNoise_method;
                    app.varNoise_manual = config.varNoise_manual;
                    app.coherence_dir = config.coherence_dir;
                    app.constellation = config.constellation;
                    app.num_looks = config.num_looks;
                    app.num_spl_method = config.num_spl_method;
                    app.spline_method = config.spline_method;
                    app.num_spl_manual = config.num_spl_manual;
                    app.lambda_method = config.lambda_method;
                    app.lambda_manual = config.lambda_manual;
                    app.coll_proc = config.coll_proc;
                    app.coll_step_est = config.coll_step_est;
                    app.minMonths = config.minMonths;
                    app.cline_resolution = config.cline_resolution;
                    app.grid_resolution = config.grid_resolution;
                    app.dtCov_STC1D = config.dtCov_STC1D;
                    app.dsCov_STC1D = config.dsCov_STC1D;
                    app.dtCov_STC2D = config.dtCov_STC2D;
                    app.dsCov_STC2D = config.dsCov_STC2D;
                    app.tCovModel_STC2D = config.tCovModel_STC2D;
                    app.tCovModel_STC1D = config.tCovModel_STC1D;
                    app.sCovModel_STC2D = config.sCovModel_STC2D;
                    app.sCovModel_STC1D = config.sCovModel_STC1D;
                    app.varNoise_DET1D = config.varNoise_DET1D;
                    app.varNoise_manual_DET1D = config.varNoise_manual_DET1D;
                    app.num_spl_method_DET1D = config.num_spl_method_DET1D;
                    app.spline_method_DET1D = config.spline_method_DET1D;
                    app.num_spl_row_manual_DET1D = config.num_spl_row_manual_DET1D;
                    app.num_spl_col_manual_DET1D = config.num_spl_col_manual_DET1D;
                    app.lambda_method_DET1D = config.lambda_method_DET1D;
                    app.lambda_manual_DET1D = config.lambda_manual_DET1D;
                    app.varNoise_DET2D = config.varNoise_DET2D;
                    app.varNoise_manual_DET2D = config.varNoise_manual_DET2D;
                    app.num_spl_method_DET2D = config.num_spl_method_DET2D;
                    app.spline_method_DET2D = config.spline_method_DET2D;
                    app.num_spl_row_manual_DET2D = config.num_spl_row_manual_DET2D;
                    app.num_spl_col_manual_DET2D = config.num_spl_col_manual_DET2D;
                    app.num_spl_t_manual_DET2D = config.num_spl_t_manual_DET2D;
                    app.lambda_method_DET2D = config.lambda_method_DET2D;
                    app.lambda_manual_DET2D = config.lambda_manual_DET2D;
                    app.useInclinedMeansST = config.useInclinedMeansST;
                    app.step_t_ST = config.step_t_ST;
                    app.polyDegreeST = config.polyDegreeST;
                    app.detrendMethod = config.detrendMethod;
                    app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
                    app.MinperiodmonthsforoutlierdetectionwithsplinesEditFieldLabel.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
                    app.TimestepEditField.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
                    app.TimestepEditFieldLabel.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
                    app.TiltedepochwiseplanesCheckBox.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
                    app.PolynomialmaxdegreeDropDown.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
                    app.PolynomialmaxdegreeDropDownLabel.Visible = ~ismember(app.procType, {'temporal', 'temporal&NNI'});
                    app.inputfilepathxlsxcsvEditField.Value = app.filepathIN;
                    app.pythoninstallationpathEditField.Value = app.pythonPath;
                    app.firstdateDatePicker.Value = app.t0IN;
                    app.removealloutputfoldersCheckBox.Value = app.flag_ouputDir;
                    app.TimestepEditField.Value = app.step_t_ST;
                    if app.useInclinedMeansST
                        app.TiltedepochwiseplanesCheckBox.Value = true;
                    else
                        app.TiltedepochwiseplanesCheckBox.Value = false;
                    end
                    app.AOIfiletypeDropDown.Value = 'Shapefile';
                    if app.flag_AOIbb
                        app.AOIfiletypeDropDown.Value = 'Bounding Box';
                    end
                    app.shapefilepathEditField.Value = app.filepathAOI;
                    app.lonminEditField.Value = app.lonMinAOI;
                    app.lonmaxEditField.Value = app.lonMaxAOI;
                    app.latminEditField.Value = app.latMinAOI;
                    app.latmaxEditField.Value = app.latMaxAOI;
                    app.projectdimensionDropDown.Value = app.projDim;
                    switch app.procType
                        case 'temporal'
                            app.processingtypeDropDown.Value = 'Temporal';
                        case 'temporal&NNI'
                            app.processingtypeDropDown.Value = 'Temporal & NNI';
                        case 'spatialDET'
                            app.processingtypeDropDown.Value = 'Spatio-temporal Deterministic';
                        case 'spatialSTC'
                            app.processingtypeDropDown.Value = 'Spatio-temporal Stochastic';
                    end
                    app.exportobsCheckBox.Value = app.flag_PSinterp;
                    app.timeseriesinterpolationatquerypointsCheckBox.Value = app.flag_tsExtr;
                    app.interpolationfilepathEditField.Value = app.filepath_EXTR;
                    app.MethodDropDown.Value = app.varNoise_method;
                    app.cohfolderEditField.Value = app.coherence_dir;
                    app.constellationDropDown.Value = app.constellation;
                    app.EditField.Value = app.markerSize;
                    app.nlooksEditField.Value = app.num_looks;
                    app.MethodDropDown_2.Value = app.num_spl_method;
                    app.IndexDropDown.Value = app.spline_method;
                    app.LambdaDropDown.Value = app.lambda_method;
                    app.MethodDropDown_3.Value = app.coll_proc;
                    app.NoiseDropDown.Value = app.varNoise_DET1D;
                    app.NoiseDropDown_2.Value = app.varNoise_DET2D;
                    app.MethodDropDown_4.Value = app.num_spl_method_DET1D;
                    app.MethodDropDown_5.Value = app.num_spl_method_DET2D;
                    app.IndexDropDown_2.Value = app.spline_method_DET1D;
                    app.IndexDropDown_3.Value = app.spline_method_DET2D;
                    app.LambdaDropDown_2.Value = app.lambda_method_DET1D;
                    app.LambdaDropDown_3.Value = app.lambda_method_DET2D;
                    if strcmp(app.detrendMethod, 'residualAtm')
                        app.ResidualatmosphericartifactsDropDown.Value = 'yes';
                    elseif strcmp(app.detrendMethod, 'cleanObs')
                        app.ResidualatmosphericartifactsDropDown.Value = 'no';
                    else
                        % Fallback for legacy files or 'no'
                        app.ResidualatmosphericartifactsDropDown.Value = 'no';
                    end
                    fields = {
                        'manualvalueEditField', 'varNoise_manual', 0; ...
                        'manualnEditField', 'num_spl_manual', 2; ...
                        'manualnEditField_2', 'lambda_manual', 0; ...
                        'manualvalueEditField_noise', 'varNoise_manual_DET1D', 0; ...
                        'manualvalueEditField_noise_2', 'varNoise_manual_DET2D', 0; ...
                        'manualnEditField_4', 'lambda_manual_DET1D', 0; ...
                        'manualnEditField_5', 'lambda_manual_DET2D', 0; ...
                        'manualvalueEditField_temporal', 'dtCov_STC1D', 0; ...
                        'manualvalueEditField_spatial', 'dsCov_STC1D', 0; ...
                        'manualvalueEditField_temporal_2', 'dtCov_STC2D', 0; ...
                        'manualvalueEditField_spatial_2', 'dsCov_STC2D', 0; ...
                        'rownEditField', 'num_spl_row_manual_DET1D', 2; ...
                        'colnEditField', 'num_spl_col_manual_DET1D', 2; ...
                        'xnEditField', 'num_spl_row_manual_DET2D', 2; ...
                        'ynEditField', 'num_spl_col_manual_DET2D', 2; ...
                        'estimationstepEditField', 'coll_step_est', 0; ...
                        'tnEditField', 'num_spl_t_manual_DET2D', 2; ...
                        'GridresolutionmEditField', 'grid_resolution', 0; ...
                        'CenterlineresolutionmEditField', 'cline_resolution', 0; ...
                        'MinperiodmonthsforoutlierdetectionwithsplinesEditField', 'minMonths', 0};
                    for i = 1:size(fields, 1)
                        propVal = app.(fields{i, 2});

                        % Check if property is scalar and not NaN (and not empty)
                        if isscalar(propVal) && ~isnan(propVal)
                            app.(fields{i, 1}).Value = propVal;
                        else
                            % Use the default placeholder (fields{i, 3}) if empty or NaN
                            app.(fields{i, 1}).Value = fields{i, 3};
                        end
                    end
                    app.ModelDropDown_temporal.Value = app.tCovModel_STC1D;
                    app.ModelDropDown_spatial.Value = app.sCovModel_STC1D;
                    app.ModelDropDown_temporal_2.Value = app.tCovModel_STC2D;
                    app.ModelDropDown_spatial_2.Value = app.sCovModel_STC2D;
                    app.PolynomialmaxdegreeDropDown.Value = app.polyDegreeST;

                    % Handle covariance fields in LoadButton as well (to fix consistency)
                    % For dtCov_STC1D
                    temp = config.dtCov_STC1D;
                    if (ischar(temp) && strcmp(temp, 'manual')) || (isnumeric(temp) && ~isnan(temp))
                        app.MethodDropDown_temporal.Value = 'manual';
                    else
                        app.MethodDropDown_temporal.Value = 'auto';
                    end
                    app.MethodDropDown_temporalValueChanged([]);
                    if isnumeric(temp) && ~isnan(temp)
                        app.manualvalueEditField_temporal.Value = temp;
                        app.manualvalueEditField_temporalValueChanged([]);
                    end

                    % For dsCov_STC1D
                    temp = config.dsCov_STC1D;
                    if (ischar(temp) && strcmp(temp, 'manual')) || (isnumeric(temp) && ~isnan(temp))
                        app.MethodDropDown_spatial.Value = 'manual';
                    else
                        app.MethodDropDown_spatial.Value = 'auto';
                    end
                    app.MethodDropDown_spatialValueChanged([]);
                    if isnumeric(temp) && ~isnan(temp)
                        app.manualvalueEditField_spatial.Value = temp;
                        app.manualvalueEditField_spatialValueChanged([]);
                    end

                    % For dtCov_STC2D
                    temp = config.dtCov_STC2D;
                    if (ischar(temp) && strcmp(temp, 'manual')) || (isnumeric(temp) && ~isnan(temp))
                        app.MethodDropDown_temporal_2.Value = 'manual';
                    else
                        app.MethodDropDown_temporal_2.Value = 'auto';
                    end
                    app.MethodDropDown_temporal_2ValueChanged([]);
                    if isnumeric(temp) && ~isnan(temp)
                        app.manualvalueEditField_temporal_2.Value = temp;
                        app.manualvalueEditField_temporal_2ValueChanged([]);
                    end

                    % For dsCov_STC2D
                    temp = config.dsCov_STC2D;
                    if (ischar(temp) && strcmp(temp, 'manual')) || (isnumeric(temp) && ~isnan(temp))
                        app.MethodDropDown_spatial_2.Value = 'manual';
                    else
                        app.MethodDropDown_spatial_2.Value = 'auto';
                    end
                    app.MethodDropDown_spatial_2ValueChanged([]);
                    if isnumeric(temp) && ~isnan(temp)
                        app.manualvalueEditField_spatial_2.Value = temp;
                        app.manualvalueEditField_spatial_2ValueChanged([]);
                    end

                    app.AOIfiletypeDropDownValueChanged([]);
                    app.timeseriesinterpolationatquerypointsCheckBoxValueChanged([]);
                    app.projectdimensionDropDownValueChanged([]);
                    app.processingtypeDropDownValueChanged([]);
                    app.MethodDropDownValueChanged([]);
                    app.MethodDropDown_2ValueChanged([]);
                    app.LambdaDropDownValueChanged([]);
                    app.MethodDropDown_3ValueChanged([]);
                    app.MethodDropDown_temporalValueChanged([]);
                    app.MethodDropDown_spatialValueChanged([]);
                    app.MethodDropDown_temporal_2ValueChanged([]);
                    app.MethodDropDown_spatial_2ValueChanged([]);
                    app.NoiseDropDownValueChanged([]);
                    app.NoiseDropDown_2ValueChanged([]);
                    app.MethodDropDown_4ValueChanged([]);
                    app.MethodDropDown_5ValueChanged([]);
                    app.LambdaDropDown_2ValueChanged([]);
                    app.LambdaDropDown_3ValueChanged([]);
                    app.ResidualatmosphericartifactsDropDownValueChanged([]);
                    % Reset output folder label
                    app.outputDir = '';
                    app.OutputfolderLabel.Text = 'Output Folder: Not set';
                    app.LampLoad.Color = 'green';
                catch
                    app.LampLoad.Color = 'red';
                end
            else
                app.LampLoad.Color = 'red';
            end
        end

        % Value changed function: MethodDropDown
        function MethodDropDownValueChanged(app, event)
            app.varNoise_method = app.MethodDropDown.Value;
            app.manualvalueEditField.Visible = strcmp(app.varNoise_method, 'manual');
            app.manualvalueEditFieldLabel.Visible = strcmp(app.varNoise_method, 'manual');
            app.cohfolderEditField.Visible = strcmp(app.varNoise_method, 'coherence');
            app.cohfolderEditFieldLabel.Visible = strcmp(app.varNoise_method, 'coherence');
            app.BrowseButton_2.Visible = strcmp(app.varNoise_method, 'coherence');
            app.constellationDropDown.Visible = strcmp(app.varNoise_method, 'coherence');
            app.constellationDropDownLabel.Visible = strcmp(app.varNoise_method, 'coherence');
            app.nlooksEditField.Visible = strcmp(app.varNoise_method, 'coherence');
            app.nlooksEditFieldLabel.Visible = strcmp(app.varNoise_method, 'coherence');
            if ~strcmp(app.varNoise_method, 'manual')
                app.varNoise_manual = NaN;
                app.manualvalueEditField.Value = 0; % Placeholder for NaN
            end
            if ~strcmp(app.varNoise_method, 'coherence')
                app.coherence_dir = '';
                app.cohfolderEditField.Value = '';
                app.num_looks = [];
                app.nlooksEditField.Value = 1; % Default to minimum valid value
            end
        end

        % Value changed function: manualvalueEditField
        function manualvalueEditFieldValueChanged(app, event)
            value = app.manualvalueEditField.Value;
            limits = app.manualvalueEditField.Limits;

            % Check if the field is empty or user intends NaN
            if isempty(value) || isnan(value) % Handle empty or NaN input
                app.varNoise_manual = NaN;
                app.manualvalueEditField.Value = 0; % Placeholder for NaN in UI
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.varNoise_manual = value;
            else
                % Revert to previous valid value or placeholder
                if isscalar(app.varNoise_manual) && ~isnan(app.varNoise_manual) && ...
                   app.varNoise_manual >= limits(1) && app.varNoise_manual <= limits(2)
                    app.manualvalueEditField.Value = app.varNoise_manual;
                else
                    app.manualvalueEditField.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Button pushed function: BrowseButton_2
        function BrowseButton_2Pushed(app, event)
            folder = uigetdir('', 'Select Coherence Directory');
            if ~isequal(folder, 0)
                app.coherence_dir = folder;
                app.cohfolderEditField.Value = folder;
            end
        end

        % Value changed function: constellationDropDown
        function constellationDropDownValueChanged(app, event)
            app.constellation = app.constellationDropDown.Value;
        end

        % Value changed function: nlooksEditField
        function nlooksEditFieldValueChanged(app, event)
            app.num_looks = app.nlooksEditField.Value;
        end

        % Value changed function: MethodDropDown_2
        function MethodDropDown_2ValueChanged(app, event)
            app.num_spl_method = app.MethodDropDown_2.Value;
            app.manualnEditField.Visible = strcmp(app.num_spl_method, 'manual');
            app.manualnEditFieldLabel.Visible = strcmp(app.num_spl_method, 'manual');
            app.IndexDropDown.Visible = strcmp(app.num_spl_method, 'auto');
            app.IndexDropDownLabel.Visible = strcmp(app.num_spl_method, 'auto');
            if ~strcmp(app.num_spl_method, 'manual')
                app.num_spl_manual = NaN;
                app.manualnEditField.Value = 2; % Placeholder for NaN (minimum spline count)
            end
        end

        % Value changed function: manualnEditField
        function manualnEditFieldValueChanged(app, event)
            value = app.manualnEditField.Value;
            limits = app.manualnEditField.Limits;

            % Check if the field is empty or user intends NaN
            if isempty(value) || isnan(value)
                app.num_spl_manual = NaN;
                app.manualnEditField.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.num_spl_manual = value;
            else
                % Revert to previous valid value or placeholder
                if isscalar(app.num_spl_manual) && ~isnan(app.num_spl_manual) && ...
                   app.num_spl_manual >= limits(1) && app.num_spl_manual <= limits(2)
                    app.manualnEditField.Value = app.num_spl_manual;
                else
                    app.manualnEditField.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: IndexDropDown
        function IndexDropDownValueChanged(app, event)
            app.spline_method = app.IndexDropDown.Value;
        end

        % Value changed function: LambdaDropDown
        function LambdaDropDownValueChanged(app, event)
            app.lambda_method = app.LambdaDropDown.Value;
            app.manualnEditField_2.Visible = strcmp(app.lambda_method, 'manual');
            app.manualnEditField_2Label.Visible = strcmp(app.lambda_method, 'manual');
            if ~strcmp(app.lambda_method, 'manual')
                app.lambda_manual = NaN;
                app.manualnEditField_2.Value = 0; % Placeholder for NaN
            end
        end

        % Value changed function: manualnEditField_2
        function manualnEditField_2ValueChanged(app, event)
            value = app.manualnEditField_2.Value;
            limits = app.manualnEditField_2.Limits;

            if isempty(value) || isnan(value)
                app.lambda_manual = NaN;
                app.manualnEditField_2.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.lambda_manual = value;
            else
                if isscalar(app.lambda_manual) && ~isnan(app.lambda_manual) && ...
                   app.lambda_manual >= limits(1) && app.lambda_manual <= limits(2)
                    app.manualnEditField_2.Value = app.lambda_manual;
                else
                    app.manualnEditField_2.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: MethodDropDown_3
        function MethodDropDown_3ValueChanged(app, event)
            app.coll_proc = app.MethodDropDown_3.Value;
            app.estimationstepEditField.Visible = strcmp(app.coll_proc, 'prediction');
            app.estimationstepEditFieldLabel.Visible = strcmp(app.coll_proc, 'prediction');
            if ~strcmp(app.coll_proc, 'prediction')
                app.coll_step_est = NaN;
                app.estimationstepEditField.Value = 0; % Placeholder for NaN
            end
        end

        % Value changed function: estimationstepEditField
        function estimationstepEditFieldValueChanged(app, event)
            value = app.estimationstepEditField.Value;
            limits = app.estimationstepEditField.Limits;

            if isempty(value) || isnan(value)
                app.coll_step_est = NaN;
                app.estimationstepEditField.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.coll_step_est = value;
            else
                if isscalar(app.coll_step_est) && ~isnan(app.coll_step_est) && ...
                   app.coll_step_est >= limits(1) && app.coll_step_est <= limits(2)
                    app.estimationstepEditField.Value = app.coll_step_est;
                else
                    app.estimationstepEditField.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function:
        % MinperiodmonthsforoutlierdetectionwithsplinesEditField
        function MinperiodmonthsforoutlierdetectionwithsplinesEditFieldValueChanged(app, event)
            value = app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Value;
            limits = app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Limits;

            if isempty(value) || isnan(value)
                app.minMonths = NaN;
                app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.minMonths = value;
            else
                if isscalar(app.minMonths) && ~isnan(app.minMonths) && ...
                   app.minMonths >= limits(1) && app.minMonths <= limits(2)
                    app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Value = app.minMonths;
                else
                    app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: CenterlineresolutionmEditField
        function CenterlineresolutionmEditFieldValueChanged(app, event)
            value = app.CenterlineresolutionmEditField.Value;
            limits = app.CenterlineresolutionmEditField.Limits;

            if isempty(value) || isnan(value)
                app.cline_resolution = NaN;
                app.CenterlineresolutionmEditField.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.cline_resolution = value;
            else
                if isscalar(app.cline_resolution) && ~isnan(app.cline_resolution) && ...
                   app.cline_resolution >= limits(1) && app.cline_resolution <= limits(2)
                    app.CenterlineresolutionmEditField.Value = app.cline_resolution;
                else
                    app.CenterlineresolutionmEditField.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: GridresolutionmEditField
        function GridresolutionmEditFieldValueChanged(app, event)
            value = app.GridresolutionmEditField.Value;
            limits = app.GridresolutionmEditField.Limits;

            if isempty(value) || isnan(value)
                app.grid_resolution = NaN;
                app.GridresolutionmEditField.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.grid_resolution = value;
            else
                if isscalar(app.grid_resolution) && ~isnan(app.grid_resolution) && ...
                   app.grid_resolution >= limits(1) && app.grid_resolution <= limits(2)
                    app.GridresolutionmEditField.Value = app.grid_resolution;
                else
                    app.GridresolutionmEditField.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: MethodDropDown_temporal
        function MethodDropDown_temporalValueChanged(app, event)
            app.dtCov_STC1D = app.MethodDropDown_temporal.Value;
            app.manualvalueEditField_temporal.Visible = strcmp(app.dtCov_STC1D, 'manual');
            app.manualvalueEditField_2Label.Visible = strcmp(app.dtCov_STC1D, 'manual');
            if ~strcmp(app.dtCov_STC1D, 'manual')
                app.dtCov_STC1D = NaN;
                app.manualvalueEditField_temporal.Value = 0; % Placeholder for NaN
            end
        end

        % Value changed function: manualvalueEditField_temporal
        function manualvalueEditField_temporalValueChanged(app, event)
            value = app.manualvalueEditField_temporal.Value;
            limits = app.manualvalueEditField_temporal.Limits;

            if isempty(value) || isnan(value)
                app.dtCov_STC1D = NaN;
                app.manualvalueEditField_temporal.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.dtCov_STC1D = value;
            else
                if isscalar(app.dtCov_STC1D) && ~isnan(app.dtCov_STC1D) && ...
                   app.dtCov_STC1D >= limits(1) && app.dtCov_STC1D <= limits(2)
                    app.manualvalueEditField_temporal.Value = app.dtCov_STC1D;
                else
                    app.manualvalueEditField_temporal.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: MethodDropDown_spatial
        function MethodDropDown_spatialValueChanged(app, event)
            app.dsCov_STC1D = app.MethodDropDown_spatial.Value;
            app.manualvalueEditField_spatial.Visible = strcmp(app.dsCov_STC1D, 'manual');
            app.manualvalueEditField_3Label.Visible = strcmp(app.dsCov_STC1D, 'manual');
            if ~strcmp(app.dsCov_STC1D, 'manual')
                app.dsCov_STC1D = NaN;
                app.manualvalueEditField_spatial.Value = 0; % Placeholder for NaN
            end
        end

        % Value changed function: manualvalueEditField_spatial
        function manualvalueEditField_spatialValueChanged(app, event)
            value = app.manualvalueEditField_spatial.Value;
            limits = app.manualvalueEditField_spatial.Limits;

            if isempty(value) || isnan(value)
                app.dsCov_STC1D = NaN;
                app.manualvalueEditField_spatial.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.dsCov_STC1D = value;
            else
                if isscalar(app.dsCov_STC1D) && ~isnan(app.dsCov_STC1D) && ...
                   app.dsCov_STC1D >= limits(1) && app.dsCov_STC1D <= limits(2)
                    app.manualvalueEditField_spatial.Value = app.dsCov_STC1D;
                else
                    app.manualvalueEditField_spatial.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: MethodDropDown_4
        function MethodDropDown_4ValueChanged(app, event)
            app.num_spl_method_DET1D = app.MethodDropDown_4.Value;
            app.rownEditField.Visible = strcmp(app.num_spl_method_DET1D, 'manual');
            app.rownEditFieldLabel.Visible = strcmp(app.num_spl_method_DET1D, 'manual');
            app.colnEditField.Visible = strcmp(app.num_spl_method_DET1D, 'manual');
            app.colnEditFieldLabel.Visible = strcmp(app.num_spl_method_DET1D, 'manual');
            app.IndexDropDown_2.Visible = strcmp(app.num_spl_method_DET1D, 'auto');
            app.IndexDropDown_2Label.Visible = strcmp(app.num_spl_method_DET1D, 'auto');
            if ~strcmp(app.num_spl_method_DET1D, 'manual')
                app.num_spl_row_manual_DET1D = NaN;
                app.num_spl_col_manual_DET1D = NaN;
                app.rownEditField.Value = 2; % Placeholder for NaN
                app.colnEditField.Value = 2; % Placeholder for NaN
            end
        end

        % Value changed function: MethodDropDown_5
        function MethodDropDown_5ValueChanged(app, event)
            app.num_spl_method_DET2D = app.MethodDropDown_5.Value;
            app.xnEditField.Visible = strcmp(app.num_spl_method_DET2D, 'manual');
            app.xnEditFieldLabel.Visible = strcmp(app.num_spl_method_DET2D, 'manual');
            app.ynEditField.Visible = strcmp(app.num_spl_method_DET2D, 'manual');
            app.ynEditFieldLabel.Visible = strcmp(app.num_spl_method_DET2D, 'manual');
            app.tnEditField.Visible = strcmp(app.num_spl_method_DET2D, 'manual');
            app.tnEditFieldLabel.Visible = strcmp(app.num_spl_method_DET2D, 'manual');
            app.IndexDropDown_3.Visible = strcmp(app.num_spl_method_DET2D, 'auto');
            app.IndexDropDown_3Label.Visible = strcmp(app.num_spl_method_DET2D, 'auto');
            if ~strcmp(app.num_spl_method_DET2D, 'manual')
                app.num_spl_row_manual_DET2D = NaN;
                app.num_spl_col_manual_DET2D = NaN;
                app.num_spl_t_manual_DET2D = NaN;
                app.xnEditField.Value = 2; % Placeholder for NaN
                app.ynEditField.Value = 2; % Placeholder for NaN
                app.tnEditField.Value = 2; % Placeholder for NaN
            end
        end

        % Value changed function: rownEditField
        function rownEditFieldValueChanged(app, event)
            value = app.rownEditField.Value;
            limits = app.rownEditField.Limits;

            if isempty(value) || isnan(value)
                app.num_spl_row_manual_DET1D = NaN;
                app.rownEditField.Value = 2; % Placeholder for NaN (minimum spline count)
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.num_spl_row_manual_DET1D = value;
            else
                if isscalar(app.num_spl_row_manual_DET1D) && ~isnan(app.num_spl_row_manual_DET1D) && ...
                   app.num_spl_row_manual_DET1D >= limits(1) && app.num_spl_row_manual_DET1D <= limits(2)
                    app.rownEditField.Value = app.num_spl_row_manual_DET1D;
                else
                    app.rownEditField.Value = 2; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: xnEditField
        function xnEditFieldValueChanged(app, event)
            value = app.xnEditField.Value;
            limits = app.xnEditField.Limits;

            if isempty(value) || isnan(value)
                app.num_spl_row_manual_DET2D = NaN;
                app.xnEditField.Value = 2; % Placeholder for NaN (minimum spline count)
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.num_spl_row_manual_DET2D = value;
            else
                if isscalar(app.num_spl_row_manual_DET2D) && ~isnan(app.num_spl_row_manual_DET2D) && ...
                   app.num_spl_row_manual_DET2D >= limits(1) && app.num_spl_row_manual_DET2D <= limits(2)
                    app.xnEditField.Value = app.num_spl_row_manual_DET2D;
                else
                    app.xnEditField.Value = 2; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: colnEditField
        function colnEditFieldValueChanged(app, event)
            value = app.colnEditField.Value;
            limits = app.colnEditField.Limits;

            if isempty(value) || isnan(value)
                app.num_spl_col_manual_DET1D = NaN;
                app.colnEditField.Value = 2; % Placeholder for NaN (minimum spline count)
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.num_spl_col_manual_DET1D = value;
            else
                if isscalar(app.num_spl_col_manual_DET1D) && ~isnan(app.num_spl_col_manual_DET1D) && ...
                   app.num_spl_col_manual_DET1D >= limits(1) && app.num_spl_col_manual_DET1D <= limits(2)
                    app.colnEditField.Value = app.num_spl_col_manual_DET1D;
                else
                    app.colnEditField.Value = 2; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: ynEditField
        function ynEditFieldValueChanged(app, event)
            value = app.ynEditField.Value;
            limits = app.ynEditField.Limits;

            if isempty(value) || isnan(value)
                app.num_spl_col_manual_DET2D = NaN;
                app.ynEditField.Value = 2; % Placeholder for NaN (minimum spline count)
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.num_spl_col_manual_DET2D = value;
            else
                if isscalar(app.num_spl_col_manual_DET2D) && ~isnan(app.num_spl_col_manual_DET2D) && ...
                   app.num_spl_col_manual_DET2D >= limits(1) && app.num_spl_col_manual_DET2D <= limits(2)
                    app.ynEditField.Value = app.num_spl_col_manual_DET2D;
                else
                    app.ynEditField.Value = 2; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: tnEditField
        function tnEditFieldValueChanged(app, event)
            value = app.tnEditField.Value;
            limits = app.tnEditField.Limits;

            if isempty(value) || isnan(value)
                app.num_spl_t_manual_DET2D = NaN;
                app.tnEditField.Value = 2; % Placeholder for NaN (minimum spline count)
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.num_spl_t_manual_DET2D = value;
            else
                if isscalar(app.num_spl_t_manual_DET2D) && ~isnan(app.num_spl_t_manual_DET2D) && ...
                   app.num_spl_t_manual_DET2D >= limits(1) && app.num_spl_t_manual_DET2D <= limits(2)
                    app.tnEditField.Value = app.num_spl_t_manual_DET2D;
                else
                    app.tnEditField.Value = 2; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: IndexDropDown_2
        function IndexDropDown_2ValueChanged(app, event)
            app.spline_method_DET1D = app.IndexDropDown_2.Value;
        end

        % Value changed function: IndexDropDown_3
        function IndexDropDown_3ValueChanged(app, event)
            app.spline_method_DET2D = app.IndexDropDown_3.Value;
        end

        % Value changed function: LambdaDropDown_2
        function LambdaDropDown_2ValueChanged(app, event)
            app.lambda_method_DET1D = app.LambdaDropDown_2.Value;
            app.manualnEditField_4.Visible = strcmp(app.lambda_method_DET1D, 'manual');
            app.manualnEditField_4Label.Visible = strcmp(app.lambda_method_DET1D, 'manual');
            if ~strcmp(app.lambda_method_DET1D, 'manual')
                app.lambda_manual_DET1D = NaN;
                app.manualnEditField_4.Value = 0; % Placeholder for NaN
            end
        end

        % Value changed function: LambdaDropDown_3
        function LambdaDropDown_3ValueChanged(app, event)
            app.lambda_method_DET2D = app.LambdaDropDown_3.Value;
            app.manualnEditField_5.Visible = strcmp(app.lambda_method_DET2D, 'manual');
            app.manualnEditField_5Label.Visible = strcmp(app.lambda_method_DET2D, 'manual');
            if ~strcmp(app.lambda_method_DET2D, 'manual')
                app.lambda_manual_DET2D = NaN;
                app.manualnEditField_5.Value = 0; % Placeholder for NaN
            end
        end

        % Value changed function: manualnEditField_4
        function manualnEditField_4ValueChanged(app, event)
            value = app.manualnEditField_4.Value;
            limits = app.manualnEditField_4.Limits;

            if isempty(value) || isnan(value)
                app.lambda_manual_DET1D = NaN;
                app.manualnEditField_4.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.lambda_manual_DET1D = value;
            else
                if isscalar(app.lambda_manual_DET1D) && ~isnan(app.lambda_manual_DET1D) && ...
                   app.lambda_manual_DET1D >= limits(1) && app.lambda_manual_DET1D <= limits(2)
                    app.manualnEditField_4.Value = app.lambda_manual_DET1D;
                else
                    app.manualnEditField_4.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: manualnEditField_5
        function manualnEditField_5ValueChanged(app, event)
            value = app.manualnEditField_5.Value;
            limits = app.manualnEditField_5.Limits;

            if isempty(value) || isnan(value)
                app.lambda_manual_DET2D = NaN;
                app.manualnEditField_5.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.lambda_manual_DET2D = value;
            else
                if isscalar(app.lambda_manual_DET2D) && ~isnan(app.lambda_manual_DET2D) && ...
                   app.lambda_manual_DET2D >= limits(1) && app.lambda_manual_DET2D <= limits(2)
                    app.manualnEditField_5.Value = app.lambda_manual_DET2D;
                else
                    app.manualnEditField_5.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: NoiseDropDown
        function NoiseDropDownValueChanged(app, event)
            app.varNoise_DET1D = app.NoiseDropDown.Value;
            app.manualvalueEditField_noise.Visible = strcmp(app.varNoise_DET1D, 'manual');
            app.manualvalueEditField_3Label_2.Visible = strcmp(app.varNoise_DET1D, 'manual');
            if ~strcmp(app.varNoise_DET1D, 'manual')
                app.varNoise_manual_DET1D = NaN;
                app.manualvalueEditField_noise.Value = 0; % Placeholder for NaN
            end
        end

        % Value changed function: NoiseDropDown_2
        function NoiseDropDown_2ValueChanged(app, event)
            app.varNoise_DET2D = app.NoiseDropDown_2.Value;
            app.manualvalueEditField_noise_2.Visible = strcmp(app.varNoise_DET2D, 'manual');
            app.manualvalueEditField_3Label_4.Visible = strcmp(app.varNoise_DET2D, 'manual');
            if ~strcmp(app.varNoise_DET2D, 'manual')
                app.varNoise_manual_DET2D = NaN;
                app.manualvalueEditField_noise_2.Value = 0; % Placeholder for NaN
            end
        end

        % Value changed function: manualvalueEditField_noise
        function manualvalueEditField_noiseValueChanged(app, event)
            value = app.manualvalueEditField_noise.Value;
            limits = app.manualvalueEditField_noise.Limits;

            if isempty(value) || isnan(value)
                app.varNoise_manual_DET1D = NaN;
                app.manualvalueEditField_noise.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.varNoise_manual_DET1D = value;
            else
                if isscalar(app.varNoise_manual_DET1D) && ~isnan(app.varNoise_manual_DET1D) && ...
                   app.varNoise_manual_DET1D >= limits(1) && app.varNoise_manual_DET1D <= limits(2)
                    app.manualvalueEditField_noise.Value = app.varNoise_manual_DET1D;
                else
                    app.manualvalueEditField_noise.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: manualvalueEditField_noise_2
        function manualvalueEditField_noise_2ValueChanged(app, event)
             value = app.manualvalueEditField_noise_2.Value;
            limits = app.manualvalueEditField_noise_2.Limits;

            if isempty(value) || isnan(value)
                app.varNoise_manual_DET2D = NaN;
                app.manualvalueEditField_noise_2.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.varNoise_manual_DET2D = value;
            else
                if isscalar(app.varNoise_manual_DET2D) && ~isnan(app.varNoise_manual_DET2D) && ...
                   app.varNoise_manual_DET2D >= limits(1) && app.varNoise_manual_DET2D <= limits(2)
                    app.manualvalueEditField_noise_2.Value = app.varNoise_manual_DET2D;
                else
                    app.manualvalueEditField_noise_2.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: MethodDropDown_spatial_2
        function MethodDropDown_spatial_2ValueChanged(app, event)
            app.dsCov_STC2D = app.MethodDropDown_spatial_2.Value;
            app.manualvalueEditField_spatial_2.Visible = strcmp(app.dsCov_STC2D, 'manual');
            app.manualvalueEditField_3Label_3.Visible = strcmp(app.dsCov_STC2D, 'manual');
            if ~strcmp(app.dsCov_STC2D, 'manual')
                app.dsCov_STC2D = NaN;
                app.manualvalueEditField_spatial_2.Value = 0; % Placeholder for NaN
            end
        end

        % Value changed function: MethodDropDown_temporal_2
        function MethodDropDown_temporal_2ValueChanged(app, event)
            app.dtCov_STC2D = app.MethodDropDown_temporal_2.Value;
            app.manualvalueEditField_temporal_2.Visible = strcmp(app.dtCov_STC2D, 'manual');
            app.manualvalueEditField_2Label_2.Visible = strcmp(app.dtCov_STC2D, 'manual');
            if ~strcmp(app.dtCov_STC2D, 'manual')
                app.dtCov_STC2D = NaN;
                app.manualvalueEditField_temporal_2.Value = 0; % Placeholder for NaN
            end
        end

        % Value changed function: manualvalueEditField_temporal_2
        function manualvalueEditField_temporal_2ValueChanged(app, event)
            value = app.manualvalueEditField_temporal_2.Value;
            limits = app.manualvalueEditField_temporal_2.Limits;

            if isempty(value) || isnan(value)
                app.dtCov_STC2D = NaN;
                app.manualvalueEditField_temporal_2.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.dtCov_STC2D = value;
            else
                if isscalar(app.dtCov_STC2D) && ~isnan(app.dtCov_STC2D) && ...
                   app.dtCov_STC2D >= limits(1) && app.dtCov_STC2D <= limits(2)
                    app.manualvalueEditField_temporal_2.Value = app.dtCov_STC2D;
                else
                    app.manualvalueEditField_temporal_2.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: manualvalueEditField_spatial_2
        function manualvalueEditField_spatial_2ValueChanged(app, event)
            value = app.manualvalueEditField_spatial_2.Value;
            limits = app.manualvalueEditField_spatial_2.Limits;

            if isempty(value) || isnan(value)
                app.dsCov_STC2D = NaN;
                app.manualvalueEditField_spatial_2.Value = 0; % Placeholder for NaN
            elseif isscalar(value) && value >= limits(1) && value <= limits(2)
                app.dsCov_STC2D = value;
            else
                if isscalar(app.dsCov_STC2D) && ~isnan(app.dsCov_STC2D) && ...
                   app.dsCov_STC2D >= limits(1) && app.dsCov_STC2D <= limits(2)
                    app.manualvalueEditField_spatial_2.Value = app.dsCov_STC2D;
                else
                    app.manualvalueEditField_spatial_2.Value = 0; % Default placeholder
                end
                uialert(app.UIFigure, sprintf('Value must be a scalar between %g and %g or empty for NaN', limits(1), limits(2)), 'Invalid Input');
            end
        end

        % Value changed function: TimestepEditField
        function TimestepEditFieldValueChanged2(app, event)
            app.step_t_ST = app.TimestepEditField.Value;
        end

        % Value changed function: ResidualatmosphericartifactsDropDown
        function ResidualatmosphericartifactsDropDownValueChanged(app, event)
            % Update internal variable based on selection
            selectedVal = app.ResidualatmosphericartifactsDropDown.Value;

            if strcmp(selectedVal, 'yes')
                app.detrendMethod = 'residualAtm';
                isCleanObs = false;
            elseif strcmp(selectedVal, 'no')
                app.detrendMethod = 'cleanObs';
                isCleanObs = true;
            end

            isSTC = strcmp(app.procType, 'spatialSTC');
            isDET = strcmp(app.procType, 'spatialDET');

            % Checkbox Visibility Logic
            % Show the checkbox if residualAtm is chosen (Applies to both STC and DET)
            app.TiltedepochwiseplanesCheckBox.Visible = ~isCleanObs && (isSTC || isDET);

            % Polynomial Degree Visibility Logic
            if isSTC
                app.PolynomialmaxdegreeDropDown.Visible = isCleanObs;
                app.PolynomialmaxdegreeDropDownLabel.Visible = isCleanObs;
            elseif isDET
                app.PolynomialmaxdegreeDropDown.Visible = true;
                app.PolynomialmaxdegreeDropDownLabel.Visible = true;
            end
        end

        % Value changed function: TiltedepochwiseplanesCheckBox
        function TiltedepochwiseplanesCheckBoxValueChanged(app, event)
            app.useInclinedMeansST = app.TiltedepochwiseplanesCheckBox.Value;
        end

        % Value changed function: ModelDropDown_temporal_2
        function ModelDropDown_temporal_2ValueChanged(app, event)
            app.tCovModel_STC2D = app.ModelDropDown_temporal_2.Value;
        end

        % Value changed function: PolynomialmaxdegreeDropDown
        function PolynomialmaxdegreeDropDownValueChanged(app, event)
            app.polyDegreeST = app.PolynomialmaxdegreeDropDown.Value;
        end

        % Value changed function: ModelDropDown_spatial_2
        function ModelDropDown_spatial_2ValueChanged(app, event)
            app.sCovModel_STC2D = app.ModelDropDown_spatial_2.Value;
        end

        % Value changed function: ModelDropDown_temporal
        function ModelDropDown_temporalValueChanged(app, event)
            app.tCovModel_STC1D = app.ModelDropDown_temporal.Value;
        end

        % Value changed function: ModelDropDown_spatial
        function ModelDropDown_spatialValueChanged(app, event)
            app.sCovModel_STC1D = app.ModelDropDown_spatial.Value;
        end

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)

            clc

            % set figures always to light mode
            set(groot, 'defaultFigureColor', 'w')     % figure background
            set(groot, 'defaultAxesColor',   'w')     % axes background
            set(groot, 'defaultAxesXColor',  'k')     % axis colors
            set(groot, 'defaultAxesYColor',  'k')
            set(groot, 'defaultAxesZColor',  'k')

            % initialize waitbar
            app.notifyBetaProgress(0, 'Preparing environment...');
            app.OutputfolderLabel.Text = 'Output Folder: Processing...';
            drawnow;

            try

                % Every legacy relative path is rooted explicitly for the
                % complete run. Results are stored beside the visible PHASE
                % shortcuts rather than inside the editable engine clone.
                runtimeRoot = phase_model_beta.projectRoot();
                previousRunFolder = pwd;
                runFolderCleanup = onCleanup(@() cd(previousRunFolder)); %#ok<NASGU>
                cd(runtimeRoot);
                outputRoot = fileparts(runtimeRoot);

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

                % Keep relative paths compatible with the scientific helpers,
                % while placing the actual result beside the PHASE shortcuts.
                outputDir = fullfile('..',outputDir);
                [created,createMessage] = mkdir(outputDir);
                if ~created
                    error('PHASE_Model_beta:outputCreateFailed', ...
                        'Could not create output folder %s: %s',outputDir,createMessage);
                end
                app.outputDir = char(java.io.File(outputDir).getCanonicalPath());
                fprintf('Output folder created: %s\n',app.outputDir);


                % folder for figures
                figsDir = fullfile(outputDir, 'figures');
                mkdir(figsDir);

                % folder for files
                filesDir = fullfile(outputDir, 'files');
                mkdir(filesDir)
                mkdir(fullfile(filesDir, 'shp'))
                mkdir(fullfile(filesDir, 'mat'))


                % - 3.3) Create processing folders for geoSplinter
                gS_input_path = fullfile(outputDir, 'geoSplinter', 'data_input');
                gS_output_path = fullfile(outputDir, 'geoSplinter', 'data_output');
                gS_job_path = fullfile(outputDir, 'geoSplinter', 'job');
                gS_synth_path = fullfile(outputDir, 'geoSplinter', 'data_synthesis');

                mkdir(gS_input_path)
                mkdir(gS_output_path)
                mkdir(gS_job_path)
                mkdir(gS_synth_path)


                % Update progress (40%)
                app.notifyBetaProgress(40, 'Preparing data...');
                drawnow;



                %% --- 4. Data preparation ---
                % Visualization of the imported data and prepration of the required
                % variables

                % - 4.1) Projected coordinates
                [xIN, yIN, utmZone] = deg2utm(lonlatIN(:,2), lonlatIN(:,1));
                utmZone = utmZone(1,:);
                xyIN = [xIN, yIN];


                % - 4.2) Filter PS inside the given AOI
                switch flag_AOIbb
                    case true
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

                    % convert coordinates to UTM
                    [xAOI, yAOI] = deg2utm(lonlatAOI(:,2), lonlatAOI(:,1));
                    xyAOI = [xAOI, yAOI];     % UTM coordinates for the bounding box

                    case false
                    % b) shapefile
                    % Use the same multipart-aware geographic geometry shown
                    % in the standalone map.
                    [lonlatAOI,~,aoiInfo] = phase_model_beta.readAoiShapefile( ...
                        filepathAOI,filepathIN);
                    fprintf('The shapefile AOI coordinates are %s.\n',aoiInfo.coordinateType);
                    fprintf('AOI polygon parts: %d\n',aoiInfo.partCount);
                    finiteAOI = all(isfinite(lonlatAOI),2);
                    xyAOI = NaN(size(lonlatAOI));
                    [xAOI,yAOI] = deg2utm( ...
                        lonlatAOI(finiteAOI,2),lonlatAOI(finiteAOI,1));
                    xyAOI(finiteAOI,:) = [xAOI,yAOI];

                end

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
                        [lat_grid_AOI, lon_grid_AOI] = utm2deg(x_grid_AOI, y_grid_AOI, repmat(utmZone, size(x_grid_AOI, 1), 1));
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
                        if any(contains(coord_columns, 'lat')) && any(contains(coord_columns, 'lon'))
                            latEXTR = coords(:, 1);    % first column (latitude)
                            lonEXTR = coords(:, 2);    % second column (longitude)
                        elseif any(contains(coord_columns, 'x')) && any(contains(coord_columns, 'y'))
                            xEXTR = coords(:, 1);      % first column (X)
                            yEXTR = coords(:, 2);      % second column (Y)
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
                        [xEXTR, yEXTR] = deg2utm(latEXTR, lonEXTR);
                        inside_AOI_flag = inpolygon(xEXTR, yEXTR, xyAOI(:,1), xyAOI(:,2));
                    elseif exist('xEXTR', 'var') && exist('yEXTR', 'var')
                        inside_AOI_flag = inpolygon(xEXTR, yEXTR, xyAOI(:,1), xyAOI(:,2));
                        [latEXTR, lonEXTR] = utm2deg(xEXTR, yEXTR, repmat(utmZone, size(yEXTR,1), 1));
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


                % - 7.6) Construct the information table for the uncertainties
                % initialize fields for the table
                avg_std = zeros(num_PS, 1);           % average velocity in mm/year

                % compute average velocity for each PS
                for i = 1:num_PS
                    p = polyfit(t_relFIN, stdAOI_FIN(i, :), 1);
                    avg_std(i) = p(1) * 365.25;
                end

                % combine all data into a table
                headers_FIN = [{'ID new', 'Longitude', 'Latitude', 'X (UTM)', 'Y (UTM)', 'Avg. std [mm]'}, ...
                    cellstr(datestr(t_dateFIN))'];

                % create the main table
                table_stdFIN = [num2cell(PSidIN_FIN), ...
                                    num2cell(lonlatIN_AOI_FIN), ...
                                    num2cell(xyIN_AOI_FIN), ...
                                    num2cell(avg_std), ...
                                    num2cell(stdAOI_FIN)];

                % add headers
                fullTable_stdFIN = [headers_FIN', table_FIN'];


                % - 7.7) Make the figure of the uncertainties over the AOI

                    % plot
                    f = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
                    geobasemap satellite;
                    hold on;
                    geoplot(lonlatAOI(:,2), lonlatAOI(:,1), 'k-', 'LineWidth', 2);
                    geoscatter(lonlatIN_AOI_FIN(:,2), lonlatIN_AOI_FIN(:,1), markerSize, avg_std, 'filled');
                    c = colorbar;
                    c.Label.String = 'Average LOS uncertainty [mm]';
                    set(gca, 'FontSize', 15)
                    title('Time-average uncertainty over AOI', 'FontSize', 20);
                    fig5_filename = fullfile(figsDir, 'AOI_uncertainty.png');
                    phase_model_beta.exportFigure(f,fig5_filename);
                    close(f);


                % - 7.8) Thresholding for warnings
                % confidence level
                alpha_sig = 5;

                % a) average velocity
                mean_avgVel = mean(avg_velocity, 'omitnan');
                std_avgVel = std(avg_velocity, 'omitnan');

                % standardize velocities and thresholds
                avgVel_standard = abs((avg_velocity - mean_avgVel) / std_avgVel);
                zlim = norminv(1 - (alpha_sig/100) / 2);
                zshift = norminv(1 - (50/100) / 2);

                % sensitivity factor for exponential curve
                a = 1.5;

                % initialize risk array
                avgVel_risk = zeros(size(avg_velocity));

                % compute risks based on standardized values
                for i = 1:length(avgVel_standard)
                    z = avgVel_standard(i);
                    if z <= zlim && z >= zshift
                        % low to moderate risk: 0% to 100%
                        avgVel_risk(i) = 100 * (1 - exp(-a * ((z - zshift) / (zlim - zshift))));
                    elseif z > zlim
                        % high risk: over 100%
                        avgVel_risk(i) = 100 + (z - zlim) * 50;
                    end
                end


                % b) cumulative displacement
                cumDisplacement = sum(abs(displAOI_FIN), 2);
                mean_cumDisplacement = mean(cumDisplacement, 'omitnan');
                std_cumDisplacement = std(cumDisplacement, 'omitnan');

                % standardize cumulative displacement
                cumD_standard = abs((cumDisplacement - mean_cumDisplacement) / std_cumDisplacement);

                % initialize risk array for cumulative displacement
                cumD_risk = zeros(size(cumDisplacement));

                % compute risks for cumulative displacement
                for i = 1:length(cumD_standard)
                    z = cumD_standard(i);
                    if z <= zlim && z >= zshift
                        % low to moderate risk: 0% to 100%
                        cumD_risk(i) = 100 * (1 - exp(-a * ((z - zshift) / (zlim - zshift))));
                    elseif z > zlim
                        % high risk: over 100%
                        cumD_risk(i) = 100 + (z - zlim) * 50;
                    end
                end


                % c) global risk
                global_risk = 0.5 * avgVel_risk + 0.5 * cumD_risk;
                risk_table = table(PSidIN_FIN, lonlatIN_AOI_FIN(:,1), lonlatIN_AOI_FIN(:,2), xyIN_AOI_FIN(:,1), xyIN_AOI_FIN(:,2), ...
                    avg_velocity, avgVel_risk, cumDisplacement, cumD_risk, global_risk, ...
                    'VariableNames', {'ID new', 'Longitude', 'Latitude', 'X (UTM)', 'Y (UTM)', 'Avg velocity [mm/yr]', ...
                    'Velocity risk [%]', 'Cumulative displ. [mm]', 'Cumulative displ. risk [%]', 'Global risk [%]'});


                % - 7.9) Shapefile creation and export
                % Check input consistency
                if size(lonlatIN_AOI_FIN, 1) ~= size(xyIN_AOI_FIN, 1) || ...
                   size(lonlatIN_AOI_FIN, 1) ~= size(PSidIN_FIN, 1) || ...
                   size(lonlatIN_AOI_FIN, 1) ~= size(avg_velocity, 1) || ...
                   size(lonlatIN_AOI_FIN, 1) ~= size(global_risk, 1) || ...
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
                                   'GlobRisk', num2cell(global_risk));        % Global risk

                % Create a FieldNameMap for .dbf field names
                FieldNameMap = containers.Map;
                FieldNameMap('PSid') = 'ID_new';
                FieldNameMap('UTM_X') = 'X_UTM';
                FieldNameMap('UTM_Y') = 'Y_UTM';
                FieldNameMap('AvgVel') = 'Avg_Vel';
                FieldNameMap('AvgStd') = 'Avg_Std';
                FieldNameMap('GlobRisk') = 'Glob_Risk';

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


                % - 7.11) Export multi-band GeoTIFF (only for 2D spatial models)
                if strcmp(projDim, '2D')
                    fprintf('Exporting multi-band GeoTIFF of the 2D spatio-temporal model...\n');

                    % 7.11.1) Dimensions
                    n_rows = size(x_grid, 1);
                    n_cols = size(x_grid, 2);
                    n_times = size(displAOI_FIN, 2);

                    % 7.11.2) Pre-allocate 3D rasters (NaN automatically acts as NoData outside the AOI)
                    tiff_matrix = NaN(n_rows, n_cols, n_times);
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

                    % 7.11.5) Extract EPSG code from utmZone string (e.g., '32 T')
                    filenameOUT_tiff_displ = fullfile(filesDir, strcat(filenameOUT, '_model.tif'));
                    filenameOUT_tiff_std   = fullfile(filesDir, strcat(filenameOUT, '_uncertainty.tif'));

                    try
                        zone_num = sscanf(utmZone, '%d');
                        if any(upper(utmZone) >= 'N')
                            epsg_code = 32600 + zone_num; % Northern hemisphere
                        else
                            epsg_code = 32700 + zone_num; % Southern hemisphere
                        end

                        geotiffwrite(filenameOUT_tiff_displ, tiff_matrix_displ, R, 'CoordRefSysCode', epsg_code);
                        geotiffwrite(filenameOUT_tiff_std, tiff_matrix_std, R, 'CoordRefSysCode', epsg_code);
                        fprintf('GeoTIFFs (Model and Uncertainty) exported successfully.\n');
                    catch ME
                        % Fallback if EPSG extraction fails
                        geotiffwrite(filenameOUT_tiff_displ, tiff_matrix_displ, R);
                        geotiffwrite(filenameOUT_tiff_std, tiff_matrix_std, R);
                        fprintf('GeoTIFFs exported (without explicit EPSG tag): %s\n', ME.message);
                    end
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
                    ['Chosen modelling: ', OUTprojDim], '';            % project spatial dimension
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
                subtitle2_s2 = 'Coordinates are given in WGS84 (geographic [deg]) and UTM (projected [m] - according to the proper fuse).';
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
                subtitle2_s3 = 'Coordinates are given in WGS84 (geographic [deg]) and UTM (projected [m] - according to the proper fuse).';
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
                subtitle1_s4 = ['Each PS displacement time series uncertainty is obtained through the ' ...
                    'contribution of all the components present in the modelled displacement.'];
                subtitle2_s4 = 'Coordinates are given in WGS84 (geographic [deg]) and UTM (projected [m] - according to the proper fuse).';
                subtitle3_s4 = 'Displacement is given along the satellite Line Of Sight (LOS) [mm]. The average velocity is in [mm/year].';
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
                sheetName5 = 'Alerts';

                % define the placeholders for the title and subtitles
                title_s5 = 'Automatic alerts';
                subtitle1_s5 = ['Alerts thresholds are computed considering both average velocity ' ...
                    'and cumulative displacement, for each PS.'];
                subtitle2_s5 = 'Threshold values are obtained considering a significance level of 5%.';
                subtitle3_s5 = 'Automatic alerts for PS - inside the AOI';

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
                writetable(risk_table, filenameOUT_e, 'Sheet', sheetName5, 'Range', startTableCell);

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
                app.OutputfolderLabel.Text = sprintf('Output Folder: %s (Processing Complete)', outputDir);
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
                app.OutputfolderLabel.Text = 'Status: Error occurred';

                % 4. hand the error to the visible standalone controller
                app.notifyBetaProgress(100, ['Processing failed: ' err.message]);
                rethrow(err);

            end

        end

        % Value changed function: lonmaxEditField
        function lonmaxEditFieldValueChanged(app, event)
            app.lonMaxAOI = app.lonmaxEditField.Value;
        end

        % Value changed function: lonminEditField
        function lonminEditFieldValueChanged(app, event)
            app.lonMinAOI = app.lonminEditField.Value;
        end

        % Value changed function: latmaxEditField
        function latmaxEditFieldValueChanged(app, event)
            app.latMaxAOI = app.latmaxEditField.Value;
        end

        % Value changed function: latminEditField
        function latminEditFieldValueChanged(app, event)
            app.latMinAOI = app.latminEditField.Value;
        end

        % Value changed function: exportobsCheckBox
        function exportobsCheckBoxValueChanged(app, event)
            app.flag_PSinterp = app.exportobsCheckBox.Value;
        end
    end

    % Component initialization
    methods (Access = public)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = phase_model_beta.projectRoot();

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.Position = [100 100 1000 650];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Theme = 'light';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 57 1000 478];

            % Create InputFilesTab
            app.InputFilesTab = uitab(app.TabGroup);
            app.InputFilesTab.Title = 'Input Files';
            app.InputFilesTab.BackgroundColor = [1 1 1];

            % Create inputfilepathxlsxcsvEditFieldLabel
            app.inputfilepathxlsxcsvEditFieldLabel = uilabel(app.InputFilesTab);
            app.inputfilepathxlsxcsvEditFieldLabel.HorizontalAlignment = 'right';
            app.inputfilepathxlsxcsvEditFieldLabel.FontName = 'Montserrat';
            app.inputfilepathxlsxcsvEditFieldLabel.FontSize = 13;
            app.inputfilepathxlsxcsvEditFieldLabel.Position = [417 361 295 22];
            app.inputfilepathxlsxcsvEditFieldLabel.Text = 'input filepath (.xlsx / .csv)';

            % Create inputfilepathxlsxcsvEditField
            app.inputfilepathxlsxcsvEditField = uieditfield(app.InputFilesTab, 'text');
            app.inputfilepathxlsxcsvEditField.ValueChangedFcn = createCallbackFcn(app, @inputfilepathxlsxcsvEditFieldValueChanged, true);
            app.inputfilepathxlsxcsvEditField.FontName = 'Montserrat';
            app.inputfilepathxlsxcsvEditField.Position = [29 359 389 26];

            % Create BrowseButton
            app.BrowseButton = uibutton(app.InputFilesTab, 'push');
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);
            app.BrowseButton.FontName = 'Montserrat';
            app.BrowseButton.Position = [435 360 100 23];
            app.BrowseButton.Text = 'Browse';

            % Create firstdateDatePickerLabel
            app.firstdateDatePickerLabel = uilabel(app.InputFilesTab);
            app.firstdateDatePickerLabel.HorizontalAlignment = 'right';
            app.firstdateDatePickerLabel.FontName = 'Montserrat';
            app.firstdateDatePickerLabel.FontSize = 13;
            app.firstdateDatePickerLabel.Position = [189 216 65 22];
            app.firstdateDatePickerLabel.Text = 'first date';

            % Create firstdateDatePicker
            app.firstdateDatePicker = uidatepicker(app.InputFilesTab);
            app.firstdateDatePicker.ValueChangedFcn = createCallbackFcn(app, @firstdateDatePickerValueChanged, true);
            app.firstdateDatePicker.Position = [29 216 150 22];

            % Create pythoninstallationpathEditFieldLabel
            app.pythoninstallationpathEditFieldLabel = uilabel(app.InputFilesTab);
            app.pythoninstallationpathEditFieldLabel.HorizontalAlignment = 'right';
            app.pythoninstallationpathEditFieldLabel.FontName = 'Montserrat';
            app.pythoninstallationpathEditFieldLabel.FontSize = 13;
            app.pythoninstallationpathEditFieldLabel.Position = [417 285 168 22];
            app.pythoninstallationpathEditFieldLabel.Text = 'python installation path';

            % Create pythoninstallationpathEditField
            app.pythoninstallationpathEditField = uieditfield(app.InputFilesTab, 'text');
            app.pythoninstallationpathEditField.ValueChangedFcn = createCallbackFcn(app, @pythoninstallationpathEditFieldValueChanged, true);
            app.pythoninstallationpathEditField.FontName = 'Montserrat';
            app.pythoninstallationpathEditField.Position = [29 283 389 26];

            % Create removealloutputfoldersCheckBox
            app.removealloutputfoldersCheckBox = uicheckbox(app.InputFilesTab);
            app.removealloutputfoldersCheckBox.ValueChangedFcn = createCallbackFcn(app, @removealloutputfoldersCheckBoxValueChanged, true);
            app.removealloutputfoldersCheckBox.Text = 'remove all output folders';
            app.removealloutputfoldersCheckBox.FontName = 'Montserrat';
            app.removealloutputfoldersCheckBox.FontSize = 13;
            app.removealloutputfoldersCheckBox.Position = [29 154 195 22];

            % Create AOISelectionTab
            app.AOISelectionTab = uitab(app.TabGroup);
            app.AOISelectionTab.Title = 'AOI Selection';
            app.AOISelectionTab.BackgroundColor = [1 1 1];

            % Create shapefilepathEditFieldLabel
            app.shapefilepathEditFieldLabel = uilabel(app.AOISelectionTab);
            app.shapefilepathEditFieldLabel.HorizontalAlignment = 'right';
            app.shapefilepathEditFieldLabel.FontName = 'Montserrat';
            app.shapefilepathEditFieldLabel.FontSize = 13;
            app.shapefilepathEditFieldLabel.Position = [417 285 229 22];
            app.shapefilepathEditFieldLabel.Text = 'shapefile path';

            % Create shapefilepathEditField
            app.shapefilepathEditField = uieditfield(app.AOISelectionTab, 'text');
            app.shapefilepathEditField.FontName = 'Montserrat';
            app.shapefilepathEditField.Position = [29 283 389 26];

            % Create AOIfiletypeDropDownLabel
            app.AOIfiletypeDropDownLabel = uilabel(app.AOISelectionTab);
            app.AOIfiletypeDropDownLabel.HorizontalAlignment = 'right';
            app.AOIfiletypeDropDownLabel.FontName = 'Montserrat';
            app.AOIfiletypeDropDownLabel.FontSize = 13;
            app.AOIfiletypeDropDownLabel.Position = [125 361 99 22];
            app.AOIfiletypeDropDownLabel.Text = 'AOI file type';

            % Create AOIfiletypeDropDown
            app.AOIfiletypeDropDown = uidropdown(app.AOISelectionTab);
            app.AOIfiletypeDropDown.Items = {'Shapefile', 'Bounding Box'};
            app.AOIfiletypeDropDown.ValueChangedFcn = createCallbackFcn(app, @AOIfiletypeDropDownValueChanged, true);
            app.AOIfiletypeDropDown.FontName = 'Montserrat';
            app.AOIfiletypeDropDown.FontSize = 13;
            app.AOIfiletypeDropDown.Position = [29 361 100 22];
            app.AOIfiletypeDropDown.Value = 'Shapefile';

            % Create latminEditFieldLabel
            app.latminEditFieldLabel = uilabel(app.AOISelectionTab);
            app.latminEditFieldLabel.HorizontalAlignment = 'right';
            app.latminEditFieldLabel.FontName = 'Montserrat';
            app.latminEditFieldLabel.FontSize = 13;
            app.latminEditFieldLabel.Position = [139 228 51 22];
            app.latminEditFieldLabel.Text = 'lat min';

            % Create latminEditField
            app.latminEditField = uieditfield(app.AOISelectionTab, 'numeric');
            app.latminEditField.Limits = [-90 90];
            app.latminEditField.AllowEmpty = 'on';
            app.latminEditField.ValueChangedFcn = createCallbackFcn(app, @latminEditFieldValueChanged, true);
            app.latminEditField.FontName = 'Montserrat';
            app.latminEditField.FontSize = 13;
            app.latminEditField.Position = [29 228 100 22];

            % Create latmaxEditFieldLabel
            app.latmaxEditFieldLabel = uilabel(app.AOISelectionTab);
            app.latmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.latmaxEditFieldLabel.FontName = 'Montserrat';
            app.latmaxEditFieldLabel.FontSize = 13;
            app.latmaxEditFieldLabel.Position = [472 228 53 22];
            app.latmaxEditFieldLabel.Text = 'lat max';

            % Create latmaxEditField
            app.latmaxEditField = uieditfield(app.AOISelectionTab, 'numeric');
            app.latmaxEditField.Limits = [-90 90];
            app.latmaxEditField.AllowEmpty = 'on';
            app.latmaxEditField.ValueChangedFcn = createCallbackFcn(app, @latmaxEditFieldValueChanged, true);
            app.latmaxEditField.FontName = 'Montserrat';
            app.latmaxEditField.FontSize = 13;
            app.latmaxEditField.Position = [364 228 100 22];

            % Create lonminEditFieldLabel
            app.lonminEditFieldLabel = uilabel(app.AOISelectionTab);
            app.lonminEditFieldLabel.HorizontalAlignment = 'right';
            app.lonminEditFieldLabel.FontName = 'Montserrat';
            app.lonminEditFieldLabel.FontSize = 13;
            app.lonminEditFieldLabel.Position = [136 285 55 22];
            app.lonminEditFieldLabel.Text = 'lon min';

            % Create lonminEditField
            app.lonminEditField = uieditfield(app.AOISelectionTab, 'numeric');
            app.lonminEditField.Limits = [-180 180];
            app.lonminEditField.AllowEmpty = 'on';
            app.lonminEditField.ValueChangedFcn = createCallbackFcn(app, @lonminEditFieldValueChanged, true);
            app.lonminEditField.FontName = 'Montserrat';
            app.lonminEditField.FontSize = 13;
            app.lonminEditField.Position = [30 285 100 22];

            % Create lonmaxEditFieldLabel
            app.lonmaxEditFieldLabel = uilabel(app.AOISelectionTab);
            app.lonmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.lonmaxEditFieldLabel.FontName = 'Montserrat';
            app.lonmaxEditFieldLabel.FontSize = 13;
            app.lonmaxEditFieldLabel.Position = [468 285 57 22];
            app.lonmaxEditFieldLabel.Text = 'lon max';

            % Create lonmaxEditField
            app.lonmaxEditField = uieditfield(app.AOISelectionTab, 'numeric');
            app.lonmaxEditField.Limits = [-180 180];
            app.lonmaxEditField.AllowEmpty = 'on';
            app.lonmaxEditField.ValueChangedFcn = createCallbackFcn(app, @lonmaxEditFieldValueChanged, true);
            app.lonmaxEditField.FontName = 'Montserrat';
            app.lonmaxEditField.FontSize = 13;
            app.lonmaxEditField.Position = [364 285 100 22];

            % Create BrowseShapefileButton
            app.BrowseShapefileButton = uibutton(app.AOISelectionTab, 'push');
            app.BrowseShapefileButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseShapefileButtonPushed, true);
            app.BrowseShapefileButton.FontName = 'Montserrat';
            app.BrowseShapefileButton.FontSize = 13;
            app.BrowseShapefileButton.Position = [437 283 100 24];
            app.BrowseShapefileButton.Text = 'Browse';

            % Create ProcessingOptionsTab
            app.ProcessingOptionsTab = uitab(app.TabGroup);
            app.ProcessingOptionsTab.Title = 'Processing Options';
            app.ProcessingOptionsTab.BackgroundColor = [1 1 1];

            % Create projectdimensionDropDownLabel
            app.projectdimensionDropDownLabel = uilabel(app.ProcessingOptionsTab);
            app.projectdimensionDropDownLabel.HorizontalAlignment = 'right';
            app.projectdimensionDropDownLabel.FontName = 'Montserrat';
            app.projectdimensionDropDownLabel.FontSize = 13;
            app.projectdimensionDropDownLabel.Position = [100 361 165 22];
            app.projectdimensionDropDownLabel.Text = 'project dimension';

            % Create projectdimensionDropDown
            app.projectdimensionDropDown = uidropdown(app.ProcessingOptionsTab);
            app.projectdimensionDropDown.Items = {'1D', '2D'};
            app.projectdimensionDropDown.ValueChangedFcn = createCallbackFcn(app, @projectdimensionDropDownValueChanged, true);
            app.projectdimensionDropDown.FontName = 'Montserrat';
            app.projectdimensionDropDown.FontSize = 13;
            app.projectdimensionDropDown.Position = [29 361 100 22];
            app.projectdimensionDropDown.Value = '1D';

            % Create processingtypeDropDownLabel
            app.processingtypeDropDownLabel = uilabel(app.ProcessingOptionsTab);
            app.processingtypeDropDownLabel.HorizontalAlignment = 'right';
            app.processingtypeDropDownLabel.FontName = 'Montserrat';
            app.processingtypeDropDownLabel.FontSize = 13;
            app.processingtypeDropDownLabel.Position = [238 293 167 22];
            app.processingtypeDropDownLabel.Text = 'processing type';

            % Create processingtypeDropDown
            app.processingtypeDropDown = uidropdown(app.ProcessingOptionsTab);
            app.processingtypeDropDown.Items = {'Temporal', 'Temporal & NNI', 'Spatio-temporal Deterministic', 'Spatio-temporal Stochastic'};
            app.processingtypeDropDown.ValueChangedFcn = createCallbackFcn(app, @processingtypeDropDownValueChanged, true);
            app.processingtypeDropDown.FontName = 'Montserrat';
            app.processingtypeDropDown.FontSize = 13;
            app.processingtypeDropDown.Position = [29 293 256 22];
            app.processingtypeDropDown.Value = 'Temporal';

            % Create EditField
            app.EditField = uieditfield(app.ProcessingOptionsTab, 'numeric');
            app.EditField.Limits = [0 Inf];
            app.EditField.RoundFractionalValues = 'on';
            app.EditField.ValueDisplayFormat = '%.0f';
            app.EditField.ValueChangedFcn = createCallbackFcn(app, @EditFieldValueChanged, true);
            app.EditField.HorizontalAlignment = 'left';
            app.EditField.FontName = 'Montserrat';
            app.EditField.FontSize = 13;
            app.EditField.Position = [30 226 99 22];
            app.EditField.Value = 10;

            % Create markersizeforplotsLabel
            app.markersizeforplotsLabel = uilabel(app.ProcessingOptionsTab);
            app.markersizeforplotsLabel.FontName = 'Montserrat';
            app.markersizeforplotsLabel.FontSize = 13;
            app.markersizeforplotsLabel.Position = [146 226 136 23];
            app.markersizeforplotsLabel.Text = 'marker size for plots';

            % Create AdvancedParametersTab
            app.AdvancedParametersTab = uitab(app.TabGroup);
            app.AdvancedParametersTab.Title = 'Advanced Parameters';
            app.AdvancedParametersTab.BackgroundColor = [1 1 1];

            % Create NoiseVariancePanel
            app.NoiseVariancePanel = uipanel(app.AdvancedParametersTab);
            app.NoiseVariancePanel.BorderColor = [0.502 0.502 0.502];
            app.NoiseVariancePanel.HighlightColor = [0.502 0.502 0.502];
            app.NoiseVariancePanel.Title = 'Noise Variance';
            app.NoiseVariancePanel.BackgroundColor = [1 1 1];
            app.NoiseVariancePanel.FontName = 'Montserrat';
            app.NoiseVariancePanel.FontWeight = 'bold';
            app.NoiseVariancePanel.FontSize = 13;
            app.NoiseVariancePanel.Position = [30 249 460 180];

            % Create MethodDropDownLabel
            app.MethodDropDownLabel = uilabel(app.NoiseVariancePanel);
            app.MethodDropDownLabel.HorizontalAlignment = 'right';
            app.MethodDropDownLabel.FontName = 'Montserrat';
            app.MethodDropDownLabel.Position = [18 123 52 22];
            app.MethodDropDownLabel.Text = 'Method';

            % Create MethodDropDown
            app.MethodDropDown = uidropdown(app.NoiseVariancePanel);
            app.MethodDropDown.Items = {'auto', 'manual', 'coherence'};
            app.MethodDropDown.ValueChangedFcn = createCallbackFcn(app, @MethodDropDownValueChanged, true);
            app.MethodDropDown.FontName = 'Montserrat';
            app.MethodDropDown.Position = [85 123 100 22];
            app.MethodDropDown.Value = 'auto';

            % Create manualvalueEditFieldLabel
            app.manualvalueEditFieldLabel = uilabel(app.NoiseVariancePanel);
            app.manualvalueEditFieldLabel.HorizontalAlignment = 'right';
            app.manualvalueEditFieldLabel.FontName = 'Montserrat';
            app.manualvalueEditFieldLabel.Position = [225 123 86 22];
            app.manualvalueEditFieldLabel.Text = 'manual value';

            % Create manualvalueEditField
            app.manualvalueEditField = uieditfield(app.NoiseVariancePanel, 'numeric');
            app.manualvalueEditField.Limits = [0 Inf];
            app.manualvalueEditField.ValueChangedFcn = createCallbackFcn(app, @manualvalueEditFieldValueChanged, true);
            app.manualvalueEditField.FontName = 'Montserrat';
            app.manualvalueEditField.Position = [326 123 100 22];

            % Create cohfolderEditFieldLabel
            app.cohfolderEditFieldLabel = uilabel(app.NoiseVariancePanel);
            app.cohfolderEditFieldLabel.HorizontalAlignment = 'right';
            app.cohfolderEditFieldLabel.FontName = 'Montserrat';
            app.cohfolderEditFieldLabel.Position = [226 123 65 22];
            app.cohfolderEditFieldLabel.Text = 'coh folder';

            % Create cohfolderEditField
            app.cohfolderEditField = uieditfield(app.NoiseVariancePanel, 'text');
            app.cohfolderEditField.FontName = 'Montserrat';
            app.cohfolderEditField.Position = [306 123 69 22];

            % Create BrowseButton_2
            app.BrowseButton_2 = uibutton(app.NoiseVariancePanel, 'push');
            app.BrowseButton_2.ButtonPushedFcn = createCallbackFcn(app, @BrowseButton_2Pushed, true);
            app.BrowseButton_2.FontName = 'Montserrat';
            app.BrowseButton_2.FontSize = 11;
            app.BrowseButton_2.Position = [388 123 56 22];
            app.BrowseButton_2.Text = 'Browse';

            % Create constellationDropDownLabel
            app.constellationDropDownLabel = uilabel(app.NoiseVariancePanel);
            app.constellationDropDownLabel.HorizontalAlignment = 'right';
            app.constellationDropDownLabel.FontName = 'Montserrat';
            app.constellationDropDownLabel.Position = [226 79 82 22];
            app.constellationDropDownLabel.Text = 'constellation';

            % Create constellationDropDown
            app.constellationDropDown = uidropdown(app.NoiseVariancePanel);
            app.constellationDropDown.Items = {'Sentinel-1', 'COSMO-SkyMed'};
            app.constellationDropDown.ValueChangedFcn = createCallbackFcn(app, @constellationDropDownValueChanged, true);
            app.constellationDropDown.FontName = 'Montserrat';
            app.constellationDropDown.Position = [339 79 100 22];
            app.constellationDropDown.Value = 'Sentinel-1';

            % Create nlooksEditFieldLabel
            app.nlooksEditFieldLabel = uilabel(app.NoiseVariancePanel);
            app.nlooksEditFieldLabel.HorizontalAlignment = 'right';
            app.nlooksEditFieldLabel.FontName = 'Montserrat';
            app.nlooksEditFieldLabel.Position = [229 35 52 22];
            app.nlooksEditFieldLabel.Text = 'n° looks';

            % Create nlooksEditField
            app.nlooksEditField = uieditfield(app.NoiseVariancePanel, 'numeric');
            app.nlooksEditField.Limits = [1 Inf];
            app.nlooksEditField.AllowEmpty = 'on';
            app.nlooksEditField.ValueChangedFcn = createCallbackFcn(app, @nlooksEditFieldValueChanged, true);
            app.nlooksEditField.FontName = 'Montserrat';
            app.nlooksEditField.Position = [339 35 100 22];
            app.nlooksEditField.Value = 1;

            % Create SplinesPanel
            app.SplinesPanel = uipanel(app.AdvancedParametersTab);
            app.SplinesPanel.BorderColor = [0.502 0.502 0.502];
            app.SplinesPanel.HighlightColor = [0.502 0.502 0.502];
            app.SplinesPanel.Title = 'Splines';
            app.SplinesPanel.BackgroundColor = [1 1 1];
            app.SplinesPanel.FontName = 'Montserrat';
            app.SplinesPanel.FontWeight = 'bold';
            app.SplinesPanel.FontSize = 13;
            app.SplinesPanel.Position = [507 249 460 180];

            % Create MethodDropDown_2Label
            app.MethodDropDown_2Label = uilabel(app.SplinesPanel);
            app.MethodDropDown_2Label.HorizontalAlignment = 'right';
            app.MethodDropDown_2Label.FontName = 'Montserrat';
            app.MethodDropDown_2Label.Position = [23 123 52 22];
            app.MethodDropDown_2Label.Text = 'Method';

            % Create MethodDropDown_2
            app.MethodDropDown_2 = uidropdown(app.SplinesPanel);
            app.MethodDropDown_2.Items = {'auto', 'manual'};
            app.MethodDropDown_2.ValueChangedFcn = createCallbackFcn(app, @MethodDropDown_2ValueChanged, true);
            app.MethodDropDown_2.FontName = 'Montserrat';
            app.MethodDropDown_2.Position = [90 123 100 22];
            app.MethodDropDown_2.Value = 'auto';

            % Create manualnEditFieldLabel
            app.manualnEditFieldLabel = uilabel(app.SplinesPanel);
            app.manualnEditFieldLabel.HorizontalAlignment = 'right';
            app.manualnEditFieldLabel.FontName = 'Montserrat';
            app.manualnEditFieldLabel.Position = [228 123 92 22];
            app.manualnEditFieldLabel.Text = 'manual n°';

            % Create manualnEditField
            app.manualnEditField = uieditfield(app.SplinesPanel, 'numeric');
            app.manualnEditField.Limits = [2 Inf];
            app.manualnEditField.ValueChangedFcn = createCallbackFcn(app, @manualnEditFieldValueChanged, true);
            app.manualnEditField.FontName = 'Montserrat';
            app.manualnEditField.Position = [329 123 100 22];
            app.manualnEditField.Value = 2;

            % Create IndexDropDownLabel
            app.IndexDropDownLabel = uilabel(app.SplinesPanel);
            app.IndexDropDownLabel.HorizontalAlignment = 'right';
            app.IndexDropDownLabel.FontName = 'Montserrat';
            app.IndexDropDownLabel.Position = [262 123 38 22];
            app.IndexDropDownLabel.Text = 'Index';

            % Create IndexDropDown
            app.IndexDropDown = uidropdown(app.SplinesPanel);
            app.IndexDropDown.Items = {'MDL', 'variance', 'F_test', 'chi2_test'};
            app.IndexDropDown.ValueChangedFcn = createCallbackFcn(app, @IndexDropDownValueChanged, true);
            app.IndexDropDown.FontName = 'Montserrat';
            app.IndexDropDown.Position = [329 123 100 22];
            app.IndexDropDown.Value = 'MDL';

            % Create LambdaDropDownLabel
            app.LambdaDropDownLabel = uilabel(app.SplinesPanel);
            app.LambdaDropDownLabel.HorizontalAlignment = 'right';
            app.LambdaDropDownLabel.FontName = 'Montserrat';
            app.LambdaDropDownLabel.Position = [23 79 55 22];
            app.LambdaDropDownLabel.Text = 'Lambda';

            % Create LambdaDropDown
            app.LambdaDropDown = uidropdown(app.SplinesPanel);
            app.LambdaDropDown.Items = {'auto', 'manual'};
            app.LambdaDropDown.ValueChangedFcn = createCallbackFcn(app, @LambdaDropDownValueChanged, true);
            app.LambdaDropDown.FontName = 'Montserrat';
            app.LambdaDropDown.Position = [90 79 100 22];
            app.LambdaDropDown.Value = 'auto';

            % Create manualnEditField_2Label
            app.manualnEditField_2Label = uilabel(app.SplinesPanel);
            app.manualnEditField_2Label.HorizontalAlignment = 'right';
            app.manualnEditField_2Label.FontName = 'Montserrat';
            app.manualnEditField_2Label.Position = [228 79 92 22];
            app.manualnEditField_2Label.Text = 'manual n°';

            % Create manualnEditField_2
            app.manualnEditField_2 = uieditfield(app.SplinesPanel, 'numeric');
            app.manualnEditField_2.Limits = [0 Inf];
            app.manualnEditField_2.ValueChangedFcn = createCallbackFcn(app, @manualnEditField_2ValueChanged, true);
            app.manualnEditField_2.FontName = 'Montserrat';
            app.manualnEditField_2.Position = [329 79 100 22];

            % Create CollocationPanel
            app.CollocationPanel = uipanel(app.AdvancedParametersTab);
            app.CollocationPanel.Title = 'Collocation';
            app.CollocationPanel.BackgroundColor = [1 1 1];
            app.CollocationPanel.FontName = 'Montserrat';
            app.CollocationPanel.FontWeight = 'bold';
            app.CollocationPanel.FontSize = 13;
            app.CollocationPanel.Position = [29 49 460 180];

            % Create MethodDropDown_3Label
            app.MethodDropDown_3Label = uilabel(app.CollocationPanel);
            app.MethodDropDown_3Label.HorizontalAlignment = 'right';
            app.MethodDropDown_3Label.FontName = 'Montserrat';
            app.MethodDropDown_3Label.Position = [19 118 52 22];
            app.MethodDropDown_3Label.Text = 'Method';

            % Create MethodDropDown_3
            app.MethodDropDown_3 = uidropdown(app.CollocationPanel);
            app.MethodDropDown_3.Items = {'filtering', 'prediction'};
            app.MethodDropDown_3.ValueChangedFcn = createCallbackFcn(app, @MethodDropDown_3ValueChanged, true);
            app.MethodDropDown_3.FontName = 'Montserrat';
            app.MethodDropDown_3.Position = [86 118 100 22];
            app.MethodDropDown_3.Value = 'filtering';

            % Create estimationstepEditFieldLabel
            app.estimationstepEditFieldLabel = uilabel(app.CollocationPanel);
            app.estimationstepEditFieldLabel.HorizontalAlignment = 'right';
            app.estimationstepEditFieldLabel.FontName = 'Montserrat';
            app.estimationstepEditFieldLabel.Position = [232 118 99 22];
            app.estimationstepEditFieldLabel.Text = 'estimation step';

            % Create estimationstepEditField
            app.estimationstepEditField = uieditfield(app.CollocationPanel, 'numeric');
            app.estimationstepEditField.Limits = [0 Inf];
            app.estimationstepEditField.ValueChangedFcn = createCallbackFcn(app, @estimationstepEditFieldValueChanged, true);
            app.estimationstepEditField.FontName = 'Montserrat';
            app.estimationstepEditField.Position = [340 118 100 22];

            % Create GeneralParametersPanel
            app.GeneralParametersPanel = uipanel(app.AdvancedParametersTab);
            app.GeneralParametersPanel.Title = 'General Parameters';
            app.GeneralParametersPanel.BackgroundColor = [1 1 1];
            app.GeneralParametersPanel.FontName = 'Montserrat';
            app.GeneralParametersPanel.FontWeight = 'bold';
            app.GeneralParametersPanel.FontSize = 13;
            app.GeneralParametersPanel.Position = [507 49 460 180];

            % Create MinperiodmonthsforoutlierdetectionwithsplinesEditFieldLabel
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditFieldLabel = uilabel(app.GeneralParametersPanel);
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditFieldLabel.HorizontalAlignment = 'right';
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditFieldLabel.FontName = 'Montserrat';
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditFieldLabel.Position = [17 88 330 22];
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditFieldLabel.Text = 'Min period (months) for outlier detection with splines';

            % Create MinperiodmonthsforoutlierdetectionwithsplinesEditField
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditField = uieditfield(app.GeneralParametersPanel, 'numeric');
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Limits = [0 Inf];
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.ValueChangedFcn = createCallbackFcn(app, @MinperiodmonthsforoutlierdetectionwithsplinesEditFieldValueChanged, true);
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.FontName = 'Montserrat';
            app.MinperiodmonthsforoutlierdetectionwithsplinesEditField.Position = [362 88 67 22];

            % Create CenterlineresolutionmEditFieldLabel
            app.CenterlineresolutionmEditFieldLabel = uilabel(app.GeneralParametersPanel);
            app.CenterlineresolutionmEditFieldLabel.HorizontalAlignment = 'right';
            app.CenterlineresolutionmEditFieldLabel.FontName = 'Montserrat';
            app.CenterlineresolutionmEditFieldLabel.Position = [23 123 154 22];
            app.CenterlineresolutionmEditFieldLabel.Text = 'Centerline resolution [m]';

            % Create CenterlineresolutionmEditField
            app.CenterlineresolutionmEditField = uieditfield(app.GeneralParametersPanel, 'numeric');
            app.CenterlineresolutionmEditField.Limits = [0 Inf];
            app.CenterlineresolutionmEditField.ValueChangedFcn = createCallbackFcn(app, @CenterlineresolutionmEditFieldValueChanged, true);
            app.CenterlineresolutionmEditField.FontName = 'Montserrat';
            app.CenterlineresolutionmEditField.Position = [362 123 67 22];

            % Create GridresolutionmEditFieldLabel
            app.GridresolutionmEditFieldLabel = uilabel(app.GeneralParametersPanel);
            app.GridresolutionmEditFieldLabel.HorizontalAlignment = 'right';
            app.GridresolutionmEditFieldLabel.FontName = 'Montserrat';
            app.GridresolutionmEditFieldLabel.Position = [23 123 117 22];
            app.GridresolutionmEditFieldLabel.Text = 'Grid resolution [m]';

            % Create GridresolutionmEditField
            app.GridresolutionmEditField = uieditfield(app.GeneralParametersPanel, 'numeric');
            app.GridresolutionmEditField.Limits = [0 Inf];
            app.GridresolutionmEditField.ValueChangedFcn = createCallbackFcn(app, @GridresolutionmEditFieldValueChanged, true);
            app.GridresolutionmEditField.FontName = 'Montserrat';
            app.GridresolutionmEditField.Position = [362 123 67 22];

            % Create ResidualatmosphericartifactsDropDownLabel
            app.ResidualatmosphericartifactsDropDownLabel = uilabel(app.GeneralParametersPanel);
            app.ResidualatmosphericartifactsDropDownLabel.HorizontalAlignment = 'right';
            app.ResidualatmosphericartifactsDropDownLabel.FontName = 'Montserrat';
            app.ResidualatmosphericartifactsDropDownLabel.Position = [26 56 184 22];
            app.ResidualatmosphericartifactsDropDownLabel.Text = 'Residual atmospheric artifacts';

            % Create ResidualatmosphericartifactsDropDown
            app.ResidualatmosphericartifactsDropDown = uidropdown(app.GeneralParametersPanel);
            app.ResidualatmosphericartifactsDropDown.Items = {'no', 'yes'};
            app.ResidualatmosphericartifactsDropDown.ValueChangedFcn = createCallbackFcn(app, @ResidualatmosphericartifactsDropDownValueChanged, true);
            app.ResidualatmosphericartifactsDropDown.FontName = 'Montserrat';
            app.ResidualatmosphericartifactsDropDown.Position = [362 56 67 22];
            app.ResidualatmosphericartifactsDropDown.Value = 'no';

            % Create PolynomialmaxdegreeDropDownLabel
            app.PolynomialmaxdegreeDropDownLabel = uilabel(app.GeneralParametersPanel);
            app.PolynomialmaxdegreeDropDownLabel.HorizontalAlignment = 'right';
            app.PolynomialmaxdegreeDropDownLabel.FontName = 'Montserrat';
            app.PolynomialmaxdegreeDropDownLabel.FontSize = 11;
            app.PolynomialmaxdegreeDropDownLabel.Position = [237 28 136 22];
            app.PolynomialmaxdegreeDropDownLabel.Text = 'Polynomial max degree';

            % Create PolynomialmaxdegreeDropDown
            app.PolynomialmaxdegreeDropDown = uidropdown(app.GeneralParametersPanel);
            app.PolynomialmaxdegreeDropDown.Items = {'1', '2', '3'};
            app.PolynomialmaxdegreeDropDown.ValueChangedFcn = createCallbackFcn(app, @PolynomialmaxdegreeDropDownValueChanged, true);
            app.PolynomialmaxdegreeDropDown.FontName = 'Montserrat';
            app.PolynomialmaxdegreeDropDown.FontSize = 11;
            app.PolynomialmaxdegreeDropDown.Position = [388 31 41 16];
            app.PolynomialmaxdegreeDropDown.Value = '1';

            % Create TiltedepochwiseplanesCheckBox
            app.TiltedepochwiseplanesCheckBox = uicheckbox(app.GeneralParametersPanel);
            app.TiltedepochwiseplanesCheckBox.ValueChangedFcn = createCallbackFcn(app, @TiltedepochwiseplanesCheckBoxValueChanged, true);
            app.TiltedepochwiseplanesCheckBox.Text = 'Tilted epoch-wise planes';
            app.TiltedepochwiseplanesCheckBox.FontName = 'Montserrat';
            app.TiltedepochwiseplanesCheckBox.FontSize = 11;
            app.TiltedepochwiseplanesCheckBox.Position = [272 28 158 22];

            % Create TimestepEditFieldLabel
            app.TimestepEditFieldLabel = uilabel(app.GeneralParametersPanel);
            app.TimestepEditFieldLabel.HorizontalAlignment = 'right';
            app.TimestepEditFieldLabel.FontName = 'Montserrat';
            app.TimestepEditFieldLabel.Position = [23 22 64 22];
            app.TimestepEditFieldLabel.Text = 'Time step';

            % Create TimestepEditField
            app.TimestepEditField = uieditfield(app.GeneralParametersPanel, 'numeric');
            app.TimestepEditField.Limits = [1 Inf];
            app.TimestepEditField.RoundFractionalValues = 'on';
            app.TimestepEditField.ValueDisplayFormat = '%.0f';
            app.TimestepEditField.ValueChangedFcn = createCallbackFcn(app, @TimestepEditFieldValueChanged2, true);
            app.TimestepEditField.FontName = 'Montserrat';
            app.TimestepEditField.Position = [115 22 49 22];
            app.TimestepEditField.Value = 1;

            % Create CovariancePanel
            app.CovariancePanel = uipanel(app.AdvancedParametersTab);
            app.CovariancePanel.Title = 'Covariance';
            app.CovariancePanel.BackgroundColor = [1 1 1];
            app.CovariancePanel.FontName = 'Montserrat';
            app.CovariancePanel.FontWeight = 'bold';
            app.CovariancePanel.FontSize = 13;
            app.CovariancePanel.Position = [12 9 460 180];

            % Create MethodDropDown_4Label
            app.MethodDropDown_4Label = uilabel(app.CovariancePanel);
            app.MethodDropDown_4Label.HorizontalAlignment = 'right';
            app.MethodDropDown_4Label.FontName = 'Montserrat';
            app.MethodDropDown_4Label.Position = [19 112 52 22];
            app.MethodDropDown_4Label.Text = 'Method';

            % Create MethodDropDown_temporal
            app.MethodDropDown_temporal = uidropdown(app.CovariancePanel);
            app.MethodDropDown_temporal.Items = {'auto', 'manual'};
            app.MethodDropDown_temporal.ValueChangedFcn = createCallbackFcn(app, @MethodDropDown_temporalValueChanged, true);
            app.MethodDropDown_temporal.FontName = 'Montserrat';
            app.MethodDropDown_temporal.Position = [86 112 100 22];
            app.MethodDropDown_temporal.Value = 'auto';

            % Create manualvalueEditField_2Label
            app.manualvalueEditField_2Label = uilabel(app.CovariancePanel);
            app.manualvalueEditField_2Label.HorizontalAlignment = 'right';
            app.manualvalueEditField_2Label.FontName = 'Montserrat';
            app.manualvalueEditField_2Label.Position = [245 112 86 22];
            app.manualvalueEditField_2Label.Text = 'manual value';

            % Create manualvalueEditField_temporal
            app.manualvalueEditField_temporal = uieditfield(app.CovariancePanel, 'numeric');
            app.manualvalueEditField_temporal.Limits = [0 Inf];
            app.manualvalueEditField_temporal.ValueChangedFcn = createCallbackFcn(app, @manualvalueEditField_temporalValueChanged, true);
            app.manualvalueEditField_temporal.FontName = 'Montserrat';
            app.manualvalueEditField_temporal.Position = [340 112 100 22];

            % Create TemporalLabel
            app.TemporalLabel = uilabel(app.CovariancePanel);
            app.TemporalLabel.FontName = 'Montserrat';
            app.TemporalLabel.FontWeight = 'bold';
            app.TemporalLabel.Position = [25 133 72 22];
            app.TemporalLabel.Text = '- Temporal';

            % Create SpatialLabel
            app.SpatialLabel = uilabel(app.CovariancePanel);
            app.SpatialLabel.FontName = 'Montserrat';
            app.SpatialLabel.FontWeight = 'bold';
            app.SpatialLabel.Position = [26 57 56 22];
            app.SpatialLabel.Text = '- Spatial';

            % Create MethodDropDown_5Label
            app.MethodDropDown_5Label = uilabel(app.CovariancePanel);
            app.MethodDropDown_5Label.HorizontalAlignment = 'right';
            app.MethodDropDown_5Label.FontName = 'Montserrat';
            app.MethodDropDown_5Label.Position = [19 36 52 22];
            app.MethodDropDown_5Label.Text = 'Method';

            % Create MethodDropDown_spatial
            app.MethodDropDown_spatial = uidropdown(app.CovariancePanel);
            app.MethodDropDown_spatial.Items = {'auto', 'manual'};
            app.MethodDropDown_spatial.ValueChangedFcn = createCallbackFcn(app, @MethodDropDown_spatialValueChanged, true);
            app.MethodDropDown_spatial.FontName = 'Montserrat';
            app.MethodDropDown_spatial.Position = [86 36 100 22];
            app.MethodDropDown_spatial.Value = 'auto';

            % Create manualvalueEditField_3Label
            app.manualvalueEditField_3Label = uilabel(app.CovariancePanel);
            app.manualvalueEditField_3Label.HorizontalAlignment = 'right';
            app.manualvalueEditField_3Label.FontName = 'Montserrat';
            app.manualvalueEditField_3Label.Position = [245 36 86 22];
            app.manualvalueEditField_3Label.Text = 'manual value';

            % Create manualvalueEditField_spatial
            app.manualvalueEditField_spatial = uieditfield(app.CovariancePanel, 'numeric');
            app.manualvalueEditField_spatial.Limits = [0 Inf];
            app.manualvalueEditField_spatial.ValueChangedFcn = createCallbackFcn(app, @manualvalueEditField_spatialValueChanged, true);
            app.manualvalueEditField_spatial.FontName = 'Montserrat';
            app.manualvalueEditField_spatial.Position = [339 36 100 22];

            % Create ModelDropDownLabel
            app.ModelDropDownLabel = uilabel(app.CovariancePanel);
            app.ModelDropDownLabel.HorizontalAlignment = 'right';
            app.ModelDropDownLabel.FontName = 'Montserrat';
            app.ModelDropDownLabel.Position = [19 83 43 22];
            app.ModelDropDownLabel.Text = 'Model';

            % Create ModelDropDown_temporal
            app.ModelDropDown_temporal = uidropdown(app.CovariancePanel);
            app.ModelDropDown_temporal.Items = {'gaussian', 'gaussian+cosine', 'exponential'};
            app.ModelDropDown_temporal.ValueChangedFcn = createCallbackFcn(app, @ModelDropDown_temporalValueChanged, true);
            app.ModelDropDown_temporal.FontName = 'Montserrat';
            app.ModelDropDown_temporal.Position = [86 83 141 22];
            app.ModelDropDown_temporal.Value = 'gaussian';

            % Create ModelDropDownLabel_2
            app.ModelDropDownLabel_2 = uilabel(app.CovariancePanel);
            app.ModelDropDownLabel_2.HorizontalAlignment = 'right';
            app.ModelDropDownLabel_2.FontName = 'Montserrat';
            app.ModelDropDownLabel_2.Position = [19 7 43 22];
            app.ModelDropDownLabel_2.Text = 'Model';

            % Create ModelDropDown_spatial
            app.ModelDropDown_spatial = uidropdown(app.CovariancePanel);
            app.ModelDropDown_spatial.Items = {'gaussian', 'gaussian+cosine', 'exponential'};
            app.ModelDropDown_spatial.ValueChangedFcn = createCallbackFcn(app, @ModelDropDown_spatialValueChanged, true);
            app.ModelDropDown_spatial.FontName = 'Montserrat';
            app.ModelDropDown_spatial.Position = [86 7 141 22];
            app.ModelDropDown_spatial.Value = 'gaussian';

            % Create SplinesPanel_spatial
            app.SplinesPanel_spatial = uipanel(app.AdvancedParametersTab);
            app.SplinesPanel_spatial.Title = 'Splines';
            app.SplinesPanel_spatial.BackgroundColor = [1 1 1];
            app.SplinesPanel_spatial.FontName = 'Montserrat';
            app.SplinesPanel_spatial.FontWeight = 'bold';
            app.SplinesPanel_spatial.FontSize = 13;
            app.SplinesPanel_spatial.Position = [450 228 460 180];

            % Create MethodDropDown_4Label_2
            app.MethodDropDown_4Label_2 = uilabel(app.SplinesPanel_spatial);
            app.MethodDropDown_4Label_2.HorizontalAlignment = 'right';
            app.MethodDropDown_4Label_2.FontName = 'Montserrat';
            app.MethodDropDown_4Label_2.Position = [23 123 52 22];
            app.MethodDropDown_4Label_2.Text = 'Method';

            % Create MethodDropDown_4
            app.MethodDropDown_4 = uidropdown(app.SplinesPanel_spatial);
            app.MethodDropDown_4.Items = {'auto', 'manual'};
            app.MethodDropDown_4.ValueChangedFcn = createCallbackFcn(app, @MethodDropDown_4ValueChanged, true);
            app.MethodDropDown_4.FontName = 'Montserrat';
            app.MethodDropDown_4.Position = [90 123 100 22];
            app.MethodDropDown_4.Value = 'auto';

            % Create rownEditFieldLabel
            app.rownEditFieldLabel = uilabel(app.SplinesPanel_spatial);
            app.rownEditFieldLabel.HorizontalAlignment = 'right';
            app.rownEditFieldLabel.FontName = 'Montserrat';
            app.rownEditFieldLabel.Position = [228 123 43 22];
            app.rownEditFieldLabel.Text = 'row n°';

            % Create rownEditField
            app.rownEditField = uieditfield(app.SplinesPanel_spatial, 'numeric');
            app.rownEditField.Limits = [2 Inf];
            app.rownEditField.ValueChangedFcn = createCallbackFcn(app, @rownEditFieldValueChanged, true);
            app.rownEditField.FontName = 'Montserrat';
            app.rownEditField.Position = [279 123 43 22];
            app.rownEditField.Value = 2;

            % Create IndexDropDown_2Label
            app.IndexDropDown_2Label = uilabel(app.SplinesPanel_spatial);
            app.IndexDropDown_2Label.HorizontalAlignment = 'right';
            app.IndexDropDown_2Label.FontName = 'Montserrat';
            app.IndexDropDown_2Label.Position = [262 123 38 22];
            app.IndexDropDown_2Label.Text = 'Index';

            % Create IndexDropDown_2
            app.IndexDropDown_2 = uidropdown(app.SplinesPanel_spatial);
            app.IndexDropDown_2.Items = {'MDL', 'variance', 'F_test', 'chi2_test'};
            app.IndexDropDown_2.ValueChangedFcn = createCallbackFcn(app, @IndexDropDown_2ValueChanged, true);
            app.IndexDropDown_2.FontName = 'Montserrat';
            app.IndexDropDown_2.Position = [329 123 100 22];
            app.IndexDropDown_2.Value = 'MDL';

            % Create LambdaDropDown_2Label
            app.LambdaDropDown_2Label = uilabel(app.SplinesPanel_spatial);
            app.LambdaDropDown_2Label.HorizontalAlignment = 'right';
            app.LambdaDropDown_2Label.FontName = 'Montserrat';
            app.LambdaDropDown_2Label.Position = [23 70 55 22];
            app.LambdaDropDown_2Label.Text = 'Lambda';

            % Create LambdaDropDown_2
            app.LambdaDropDown_2 = uidropdown(app.SplinesPanel_spatial);
            app.LambdaDropDown_2.Items = {'auto', 'manual'};
            app.LambdaDropDown_2.ValueChangedFcn = createCallbackFcn(app, @LambdaDropDown_2ValueChanged, true);
            app.LambdaDropDown_2.FontName = 'Montserrat';
            app.LambdaDropDown_2.Position = [90 70 100 22];
            app.LambdaDropDown_2.Value = 'auto';

            % Create manualnEditField_4Label
            app.manualnEditField_4Label = uilabel(app.SplinesPanel_spatial);
            app.manualnEditField_4Label.HorizontalAlignment = 'right';
            app.manualnEditField_4Label.FontName = 'Montserrat';
            app.manualnEditField_4Label.Position = [228 70 92 22];
            app.manualnEditField_4Label.Text = 'manual n°';

            % Create manualnEditField_4
            app.manualnEditField_4 = uieditfield(app.SplinesPanel_spatial, 'numeric');
            app.manualnEditField_4.Limits = [0 Inf];
            app.manualnEditField_4.ValueChangedFcn = createCallbackFcn(app, @manualnEditField_4ValueChanged, true);
            app.manualnEditField_4.FontName = 'Montserrat';
            app.manualnEditField_4.Position = [329 70 100 22];

            % Create MethodDropDown_5Label_2
            app.MethodDropDown_5Label_2 = uilabel(app.SplinesPanel_spatial);
            app.MethodDropDown_5Label_2.HorizontalAlignment = 'right';
            app.MethodDropDown_5Label_2.FontName = 'Montserrat';
            app.MethodDropDown_5Label_2.Position = [23 16 39 22];
            app.MethodDropDown_5Label_2.Text = 'Noise';

            % Create NoiseDropDown
            app.NoiseDropDown = uidropdown(app.SplinesPanel_spatial);
            app.NoiseDropDown.Items = {'auto', 'manual'};
            app.NoiseDropDown.ValueChangedFcn = createCallbackFcn(app, @NoiseDropDownValueChanged, true);
            app.NoiseDropDown.FontName = 'Montserrat';
            app.NoiseDropDown.Position = [90 16 100 22];
            app.NoiseDropDown.Value = 'auto';

            % Create manualvalueEditField_3Label_2
            app.manualvalueEditField_3Label_2 = uilabel(app.SplinesPanel_spatial);
            app.manualvalueEditField_3Label_2.HorizontalAlignment = 'right';
            app.manualvalueEditField_3Label_2.FontName = 'Montserrat';
            app.manualvalueEditField_3Label_2.Position = [235 16 86 22];
            app.manualvalueEditField_3Label_2.Text = 'manual value';

            % Create manualvalueEditField_noise
            app.manualvalueEditField_noise = uieditfield(app.SplinesPanel_spatial, 'numeric');
            app.manualvalueEditField_noise.Limits = [0 Inf];
            app.manualvalueEditField_noise.ValueChangedFcn = createCallbackFcn(app, @manualvalueEditField_noiseValueChanged, true);
            app.manualvalueEditField_noise.FontName = 'Montserrat';
            app.manualvalueEditField_noise.Position = [329 16 100 22];

            % Create colnEditFieldLabel
            app.colnEditFieldLabel = uilabel(app.SplinesPanel_spatial);
            app.colnEditFieldLabel.HorizontalAlignment = 'right';
            app.colnEditFieldLabel.FontName = 'Montserrat';
            app.colnEditFieldLabel.Position = [335 123 43 22];
            app.colnEditFieldLabel.Text = 'col n°';

            % Create colnEditField
            app.colnEditField = uieditfield(app.SplinesPanel_spatial, 'numeric');
            app.colnEditField.Limits = [2 Inf];
            app.colnEditField.ValueChangedFcn = createCallbackFcn(app, @colnEditFieldValueChanged, true);
            app.colnEditField.FontName = 'Montserrat';
            app.colnEditField.Position = [386 123 43 22];
            app.colnEditField.Value = 2;

            % Create CovariancePanel_2D
            app.CovariancePanel_2D = uipanel(app.AdvancedParametersTab);
            app.CovariancePanel_2D.Title = 'Covariance';
            app.CovariancePanel_2D.BackgroundColor = [1 1 1];
            app.CovariancePanel_2D.FontName = 'Montserrat';
            app.CovariancePanel_2D.FontWeight = 'bold';
            app.CovariancePanel_2D.FontSize = 13;
            app.CovariancePanel_2D.Position = [421 244 460 180];

            % Create MethodDropDown_4Label_3
            app.MethodDropDown_4Label_3 = uilabel(app.CovariancePanel_2D);
            app.MethodDropDown_4Label_3.HorizontalAlignment = 'right';
            app.MethodDropDown_4Label_3.FontName = 'Montserrat';
            app.MethodDropDown_4Label_3.Position = [19 111 52 22];
            app.MethodDropDown_4Label_3.Text = 'Method';

            % Create MethodDropDown_temporal_2
            app.MethodDropDown_temporal_2 = uidropdown(app.CovariancePanel_2D);
            app.MethodDropDown_temporal_2.Items = {'auto', 'manual'};
            app.MethodDropDown_temporal_2.ValueChangedFcn = createCallbackFcn(app, @MethodDropDown_temporal_2ValueChanged, true);
            app.MethodDropDown_temporal_2.FontName = 'Montserrat';
            app.MethodDropDown_temporal_2.Position = [86 111 100 22];
            app.MethodDropDown_temporal_2.Value = 'auto';

            % Create manualvalueEditField_2Label_2
            app.manualvalueEditField_2Label_2 = uilabel(app.CovariancePanel_2D);
            app.manualvalueEditField_2Label_2.HorizontalAlignment = 'right';
            app.manualvalueEditField_2Label_2.FontName = 'Montserrat';
            app.manualvalueEditField_2Label_2.Position = [245 111 86 22];
            app.manualvalueEditField_2Label_2.Text = 'manual value';

            % Create manualvalueEditField_temporal_2
            app.manualvalueEditField_temporal_2 = uieditfield(app.CovariancePanel_2D, 'numeric');
            app.manualvalueEditField_temporal_2.Limits = [0 Inf];
            app.manualvalueEditField_temporal_2.ValueChangedFcn = createCallbackFcn(app, @manualvalueEditField_temporal_2ValueChanged, true);
            app.manualvalueEditField_temporal_2.FontName = 'Montserrat';
            app.manualvalueEditField_temporal_2.Position = [340 111 100 22];

            % Create TemporalLabel_2
            app.TemporalLabel_2 = uilabel(app.CovariancePanel_2D);
            app.TemporalLabel_2.FontName = 'Montserrat';
            app.TemporalLabel_2.FontWeight = 'bold';
            app.TemporalLabel_2.Position = [25 132 72 22];
            app.TemporalLabel_2.Text = '- Temporal';

            % Create SpatialLabel_2
            app.SpatialLabel_2 = uilabel(app.CovariancePanel_2D);
            app.SpatialLabel_2.FontName = 'Montserrat';
            app.SpatialLabel_2.FontWeight = 'bold';
            app.SpatialLabel_2.Position = [26 59 56 22];
            app.SpatialLabel_2.Text = '- Spatial';

            % Create MethodDropDown_5Label_3
            app.MethodDropDown_5Label_3 = uilabel(app.CovariancePanel_2D);
            app.MethodDropDown_5Label_3.HorizontalAlignment = 'right';
            app.MethodDropDown_5Label_3.FontName = 'Montserrat';
            app.MethodDropDown_5Label_3.Position = [19 38 52 22];
            app.MethodDropDown_5Label_3.Text = 'Method';

            % Create MethodDropDown_spatial_2
            app.MethodDropDown_spatial_2 = uidropdown(app.CovariancePanel_2D);
            app.MethodDropDown_spatial_2.Items = {'auto', 'manual'};
            app.MethodDropDown_spatial_2.ValueChangedFcn = createCallbackFcn(app, @MethodDropDown_spatial_2ValueChanged, true);
            app.MethodDropDown_spatial_2.FontName = 'Montserrat';
            app.MethodDropDown_spatial_2.Position = [86 38 100 22];
            app.MethodDropDown_spatial_2.Value = 'auto';

            % Create manualvalueEditField_3Label_3
            app.manualvalueEditField_3Label_3 = uilabel(app.CovariancePanel_2D);
            app.manualvalueEditField_3Label_3.HorizontalAlignment = 'right';
            app.manualvalueEditField_3Label_3.FontName = 'Montserrat';
            app.manualvalueEditField_3Label_3.Position = [245 38 86 22];
            app.manualvalueEditField_3Label_3.Text = 'manual value';

            % Create manualvalueEditField_spatial_2
            app.manualvalueEditField_spatial_2 = uieditfield(app.CovariancePanel_2D, 'numeric');
            app.manualvalueEditField_spatial_2.Limits = [0 Inf];
            app.manualvalueEditField_spatial_2.ValueChangedFcn = createCallbackFcn(app, @manualvalueEditField_spatial_2ValueChanged, true);
            app.manualvalueEditField_spatial_2.FontName = 'Montserrat';
            app.manualvalueEditField_spatial_2.Position = [339 38 100 22];

            % Create ModelDropDownLabel_3
            app.ModelDropDownLabel_3 = uilabel(app.CovariancePanel_2D);
            app.ModelDropDownLabel_3.HorizontalAlignment = 'right';
            app.ModelDropDownLabel_3.FontName = 'Montserrat';
            app.ModelDropDownLabel_3.Position = [19 82 43 22];
            app.ModelDropDownLabel_3.Text = 'Model';

            % Create ModelDropDown_temporal_2
            app.ModelDropDown_temporal_2 = uidropdown(app.CovariancePanel_2D);
            app.ModelDropDown_temporal_2.Items = {'gaussian', 'gaussian+cosine', 'exponential'};
            app.ModelDropDown_temporal_2.ValueChangedFcn = createCallbackFcn(app, @ModelDropDown_temporal_2ValueChanged, true);
            app.ModelDropDown_temporal_2.FontName = 'Montserrat';
            app.ModelDropDown_temporal_2.Position = [86 82 141 22];
            app.ModelDropDown_temporal_2.Value = 'gaussian';

            % Create ModelDropDownLabel_4
            app.ModelDropDownLabel_4 = uilabel(app.CovariancePanel_2D);
            app.ModelDropDownLabel_4.HorizontalAlignment = 'right';
            app.ModelDropDownLabel_4.FontName = 'Montserrat';
            app.ModelDropDownLabel_4.Position = [19 8 43 22];
            app.ModelDropDownLabel_4.Text = 'Model';

            % Create ModelDropDown_spatial_2
            app.ModelDropDown_spatial_2 = uidropdown(app.CovariancePanel_2D);
            app.ModelDropDown_spatial_2.Items = {'gaussian', 'gaussian+cosine', 'exponential'};
            app.ModelDropDown_spatial_2.ValueChangedFcn = createCallbackFcn(app, @ModelDropDown_spatial_2ValueChanged, true);
            app.ModelDropDown_spatial_2.FontName = 'Montserrat';
            app.ModelDropDown_spatial_2.Position = [86 8 141 22];
            app.ModelDropDown_spatial_2.Value = 'gaussian';

            % Create SplinesPanel_spatial_2D
            app.SplinesPanel_spatial_2D = uipanel(app.AdvancedParametersTab);
            app.SplinesPanel_spatial_2D.Title = 'Splines';
            app.SplinesPanel_spatial_2D.BackgroundColor = [1 1 1];
            app.SplinesPanel_spatial_2D.FontName = 'Montserrat';
            app.SplinesPanel_spatial_2D.FontWeight = 'bold';
            app.SplinesPanel_spatial_2D.FontSize = 13;
            app.SplinesPanel_spatial_2D.Position = [-9 203 460 180];

            % Create MethodDropDown_5Label_5
            app.MethodDropDown_5Label_5 = uilabel(app.SplinesPanel_spatial_2D);
            app.MethodDropDown_5Label_5.HorizontalAlignment = 'right';
            app.MethodDropDown_5Label_5.FontName = 'Montserrat';
            app.MethodDropDown_5Label_5.Position = [23 123 52 22];
            app.MethodDropDown_5Label_5.Text = 'Method';

            % Create MethodDropDown_5
            app.MethodDropDown_5 = uidropdown(app.SplinesPanel_spatial_2D);
            app.MethodDropDown_5.Items = {'auto', 'manual'};
            app.MethodDropDown_5.ValueChangedFcn = createCallbackFcn(app, @MethodDropDown_5ValueChanged, true);
            app.MethodDropDown_5.FontName = 'Montserrat';
            app.MethodDropDown_5.Position = [90 123 100 22];
            app.MethodDropDown_5.Value = 'auto';

            % Create xnEditFieldLabel
            app.xnEditFieldLabel = uilabel(app.SplinesPanel_spatial_2D);
            app.xnEditFieldLabel.HorizontalAlignment = 'right';
            app.xnEditFieldLabel.FontName = 'Montserrat';
            app.xnEditFieldLabel.Position = [206 123 43 22];
            app.xnEditFieldLabel.Text = 'x n°';

            % Create xnEditField
            app.xnEditField = uieditfield(app.SplinesPanel_spatial_2D, 'numeric');
            app.xnEditField.Limits = [2 Inf];
            app.xnEditField.ValueChangedFcn = createCallbackFcn(app, @xnEditFieldValueChanged, true);
            app.xnEditField.FontName = 'Montserrat';
            app.xnEditField.Position = [255 123 22 22];
            app.xnEditField.Value = 2;

            % Create IndexDropDown_3Label
            app.IndexDropDown_3Label = uilabel(app.SplinesPanel_spatial_2D);
            app.IndexDropDown_3Label.HorizontalAlignment = 'right';
            app.IndexDropDown_3Label.FontName = 'Montserrat';
            app.IndexDropDown_3Label.Position = [262 123 38 22];
            app.IndexDropDown_3Label.Text = 'Index';

            % Create IndexDropDown_3
            app.IndexDropDown_3 = uidropdown(app.SplinesPanel_spatial_2D);
            app.IndexDropDown_3.Items = {'MDL', 'variance', 'F_test', 'chi2_test'};
            app.IndexDropDown_3.ValueChangedFcn = createCallbackFcn(app, @IndexDropDown_3ValueChanged, true);
            app.IndexDropDown_3.FontName = 'Montserrat';
            app.IndexDropDown_3.Position = [329 123 100 22];
            app.IndexDropDown_3.Value = 'MDL';

            % Create LambdaDropDown_3Label
            app.LambdaDropDown_3Label = uilabel(app.SplinesPanel_spatial_2D);
            app.LambdaDropDown_3Label.HorizontalAlignment = 'right';
            app.LambdaDropDown_3Label.FontName = 'Montserrat';
            app.LambdaDropDown_3Label.Position = [23 71 55 22];
            app.LambdaDropDown_3Label.Text = 'Lambda';

            % Create LambdaDropDown_3
            app.LambdaDropDown_3 = uidropdown(app.SplinesPanel_spatial_2D);
            app.LambdaDropDown_3.Items = {'auto', 'manual'};
            app.LambdaDropDown_3.ValueChangedFcn = createCallbackFcn(app, @LambdaDropDown_3ValueChanged, true);
            app.LambdaDropDown_3.FontName = 'Montserrat';
            app.LambdaDropDown_3.Position = [90 71 100 22];
            app.LambdaDropDown_3.Value = 'auto';

            % Create manualnEditField_5Label
            app.manualnEditField_5Label = uilabel(app.SplinesPanel_spatial_2D);
            app.manualnEditField_5Label.HorizontalAlignment = 'right';
            app.manualnEditField_5Label.FontName = 'Montserrat';
            app.manualnEditField_5Label.Position = [228 71 92 22];
            app.manualnEditField_5Label.Text = 'manual n°';

            % Create manualnEditField_5
            app.manualnEditField_5 = uieditfield(app.SplinesPanel_spatial_2D, 'numeric');
            app.manualnEditField_5.Limits = [0 Inf];
            app.manualnEditField_5.ValueChangedFcn = createCallbackFcn(app, @manualnEditField_5ValueChanged, true);
            app.manualnEditField_5.FontName = 'Montserrat';
            app.manualnEditField_5.Position = [329 71 100 22];

            % Create NoiseDropDown_2Label
            app.NoiseDropDown_2Label = uilabel(app.SplinesPanel_spatial_2D);
            app.NoiseDropDown_2Label.HorizontalAlignment = 'right';
            app.NoiseDropDown_2Label.FontName = 'Montserrat';
            app.NoiseDropDown_2Label.Position = [23 20 39 22];
            app.NoiseDropDown_2Label.Text = 'Noise';

            % Create NoiseDropDown_2
            app.NoiseDropDown_2 = uidropdown(app.SplinesPanel_spatial_2D);
            app.NoiseDropDown_2.Items = {'auto', 'manual'};
            app.NoiseDropDown_2.ValueChangedFcn = createCallbackFcn(app, @NoiseDropDown_2ValueChanged, true);
            app.NoiseDropDown_2.FontName = 'Montserrat';
            app.NoiseDropDown_2.Position = [90 20 100 22];
            app.NoiseDropDown_2.Value = 'auto';

            % Create manualvalueEditField_3Label_4
            app.manualvalueEditField_3Label_4 = uilabel(app.SplinesPanel_spatial_2D);
            app.manualvalueEditField_3Label_4.HorizontalAlignment = 'right';
            app.manualvalueEditField_3Label_4.FontName = 'Montserrat';
            app.manualvalueEditField_3Label_4.Position = [235 20 86 22];
            app.manualvalueEditField_3Label_4.Text = 'manual value';

            % Create manualvalueEditField_noise_2
            app.manualvalueEditField_noise_2 = uieditfield(app.SplinesPanel_spatial_2D, 'numeric');
            app.manualvalueEditField_noise_2.Limits = [0 Inf];
            app.manualvalueEditField_noise_2.ValueChangedFcn = createCallbackFcn(app, @manualvalueEditField_noise_2ValueChanged, true);
            app.manualvalueEditField_noise_2.FontName = 'Montserrat';
            app.manualvalueEditField_noise_2.Position = [329 20 100 22];

            % Create ynEditFieldLabel
            app.ynEditFieldLabel = uilabel(app.SplinesPanel_spatial_2D);
            app.ynEditFieldLabel.HorizontalAlignment = 'right';
            app.ynEditFieldLabel.FontName = 'Montserrat';
            app.ynEditFieldLabel.Position = [288 123 43 22];
            app.ynEditFieldLabel.Text = 'y n°';

            % Create ynEditField
            app.ynEditField = uieditfield(app.SplinesPanel_spatial_2D, 'numeric');
            app.ynEditField.Limits = [2 Inf];
            app.ynEditField.ValueChangedFcn = createCallbackFcn(app, @ynEditFieldValueChanged, true);
            app.ynEditField.FontName = 'Montserrat';
            app.ynEditField.Position = [337 123 22 22];
            app.ynEditField.Value = 2;

            % Create tnEditFieldLabel
            app.tnEditFieldLabel = uilabel(app.SplinesPanel_spatial_2D);
            app.tnEditFieldLabel.HorizontalAlignment = 'right';
            app.tnEditFieldLabel.FontName = 'Montserrat';
            app.tnEditFieldLabel.Position = [358 123 43 22];
            app.tnEditFieldLabel.Text = 't n°';

            % Create tnEditField
            app.tnEditField = uieditfield(app.SplinesPanel_spatial_2D, 'numeric');
            app.tnEditField.Limits = [2 Inf];
            app.tnEditField.ValueChangedFcn = createCallbackFcn(app, @tnEditFieldValueChanged, true);
            app.tnEditField.FontName = 'Montserrat';
            app.tnEditField.Position = [407 123 22 22];
            app.tnEditField.Value = 2;

            % Create QuerypointsTab
            app.QuerypointsTab = uitab(app.TabGroup);
            app.QuerypointsTab.Title = 'Query points';
            app.QuerypointsTab.BackgroundColor = [1 1 1];

            % Create interpolationfilepathEditFieldLabel
            app.interpolationfilepathEditFieldLabel = uilabel(app.QuerypointsTab);
            app.interpolationfilepathEditFieldLabel.HorizontalAlignment = 'right';
            app.interpolationfilepathEditFieldLabel.FontName = 'Montserrat';
            app.interpolationfilepathEditFieldLabel.FontSize = 13;
            app.interpolationfilepathEditFieldLabel.Position = [535 244 156 22];
            app.interpolationfilepathEditFieldLabel.Text = 'interpolation filepath';

            % Create interpolationfilepathEditField
            app.interpolationfilepathEditField = uieditfield(app.QuerypointsTab, 'text');
            app.interpolationfilepathEditField.FontName = 'Montserrat';
            app.interpolationfilepathEditField.Position = [29 242 389 26];

            % Create timeseriesinterpolationatquerypointsCheckBox
            app.timeseriesinterpolationatquerypointsCheckBox = uicheckbox(app.QuerypointsTab);
            app.timeseriesinterpolationatquerypointsCheckBox.ValueChangedFcn = createCallbackFcn(app, @timeseriesinterpolationatquerypointsCheckBoxValueChanged, true);
            app.timeseriesinterpolationatquerypointsCheckBox.Text = 'time series interpolation at query points';
            app.timeseriesinterpolationatquerypointsCheckBox.FontName = 'Montserrat';
            app.timeseriesinterpolationatquerypointsCheckBox.FontSize = 13;
            app.timeseriesinterpolationatquerypointsCheckBox.Position = [29 319 285 22];

            % Create BrowseExtrapolationButton
            app.BrowseExtrapolationButton = uibutton(app.QuerypointsTab, 'push');
            app.BrowseExtrapolationButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseExtrapolationButtonPushed, true);
            app.BrowseExtrapolationButton.FontName = 'Montserrat';
            app.BrowseExtrapolationButton.FontSize = 13;
            app.BrowseExtrapolationButton.Position = [437 242 100 24];
            app.BrowseExtrapolationButton.Text = 'Browse';

            % Create exportobsCheckBox
            app.exportobsCheckBox = uicheckbox(app.QuerypointsTab);
            app.exportobsCheckBox.ValueChangedFcn = createCallbackFcn(app, @exportobsCheckBoxValueChanged, true);
            app.exportobsCheckBox.Text = 'time series interpolation at input PS';
            app.exportobsCheckBox.FontName = 'Montserrat';
            app.exportobsCheckBox.FontSize = 13;
            app.exportobsCheckBox.Position = [29 360 259 22];

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Position = [13 548 225 92];
            app.Image.ImageSource = fullfile(pathToMLAPP, 'PHASE_mod2.png');

            % Create Image_2
            app.Image_2 = uiimage(app.UIFigure);
            app.Image_2.Position = [791 548 192 92];
            app.Image_2.ImageSource = fullfile(pathToMLAPP, 'PHASE_logo.png');

            % Create StartButton
            app.StartButton = uibutton(app.UIFigure, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.BackgroundColor = [0.3098 0.4353 0.8235];
            app.StartButton.FontName = 'Montserrat';
            app.StartButton.FontSize = 18;
            app.StartButton.FontWeight = 'bold';
            app.StartButton.FontColor = [1 1 1];
            app.StartButton.Position = [451 579 100 31];
            app.StartButton.Text = 'Start';

            % Create SaveButton
            app.SaveButton = uibutton(app.UIFigure, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.SaveButton.FontName = 'Montserrat';
            app.SaveButton.FontSize = 14;
            app.SaveButton.FontWeight = 'bold';
            app.SaveButton.Position = [62 15 100 25];
            app.SaveButton.Text = 'Save';

            % Create LoadButton
            app.LoadButton = uibutton(app.UIFigure, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.LoadButton.FontName = 'Montserrat';
            app.LoadButton.FontSize = 14;
            app.LoadButton.FontWeight = 'bold';
            app.LoadButton.Position = [840 15 100 25];
            app.LoadButton.Text = 'Load';

            % Create LampSave
            app.LampSave = uilamp(app.UIFigure);
            app.LampSave.Position = [30 17 20 20];

            % Create LampLoad
            app.LampLoad = uilamp(app.UIFigure);
            app.LampLoad.Position = [951 17 20 20];

            % Create OutputfolderLabel
            app.OutputfolderLabel = uilabel(app.UIFigure);
            app.OutputfolderLabel.HorizontalAlignment = 'center';
            app.OutputfolderLabel.FontName = 'Montserrat';
            app.OutputfolderLabel.FontSize = 14;
            app.OutputfolderLabel.Position = [168 15 666 24];
            app.OutputfolderLabel.Text = 'Output folder: Not set';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'off';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = LegacyEngine

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
