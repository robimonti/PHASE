% GENERATED FROM PHASE_Preprocessing.mlapp.
% Run tools/extract_preprocessing_beta_engine.py after changing the stable app.
% The complete formerly embedded App Designer source follows as editable text.

classdef LegacyEngine < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        Image2                          matlab.ui.control.Image
        Image                           matlab.ui.control.Image
        Label                           matlab.ui.control.Label
        ConstellationSwitch             matlab.ui.control.Switch
        Sentinel1Panel                  matlab.ui.container.Panel
        TabGroup                        matlab.ui.container.TabGroup
        DownloadTab                     matlab.ui.container.Tab
        PreviewPanel                    matlab.ui.container.Panel
        DownloadProgressLabel           matlab.ui.control.Label
        DownloadProgressGauge           matlab.ui.control.LinearGauge
        RecommendedFramePathLabel       matlab.ui.control.Label
        ToggleFootprintsButton          matlab.ui.control.Button
        DownloadAllButton               matlab.ui.control.Button
        SelectAllCheckBox               matlab.ui.control.CheckBox
        PreviewSelectedLabel            matlab.ui.control.Label
        DownloadSelectedButton          matlab.ui.control.Button
        PreviewResultsTable             matlab.ui.control.Table
        DownloaderToolbarPanel          matlab.ui.container.Panel
        TogglePreviewButton             matlab.ui.control.Button
        LoadLastDownloadButton          matlab.ui.control.Button
        SearchASFButton                 matlab.ui.control.Button
        ClearAOIButton                  matlab.ui.control.Button
        LoadShapefileButton             matlab.ui.control.Button
        DrawPolygonButton               matlab.ui.control.Button
        DrawRectangleButton             matlab.ui.control.Button
        ShowFiltersButton               matlab.ui.control.Button
        SignOutButton                   matlab.ui.control.Button
        SignedInLabel                   matlab.ui.control.Label
        LoginFeedbackLabel              matlab.ui.control.Label
        LoginButton                     matlab.ui.control.Button
        EarthdataPasswordEditField      matlab.ui.control.EditField
        EarthdataUsernameEditField      matlab.ui.control.EditField
        MaxLatitudeEditField            matlab.ui.control.NumericEditField
        MaxLatitudeEditFieldLabel       matlab.ui.control.Label
        MinLatitudeEditField            matlab.ui.control.NumericEditField
        MinLatitudeEditFieldLabel       matlab.ui.control.Label
        MaxLongitudeEditField           matlab.ui.control.NumericEditField
        MaxLongitudeEditFieldLabel      matlab.ui.control.Label
        MinLongitudeEditField           matlab.ui.control.NumericEditField
        MinLongitudeEditFieldLabel      matlab.ui.control.Label
        FilterPanel                     matlab.ui.container.Panel
        SamplingUnitDropDown            matlab.ui.control.DropDown
        SamplingRateEditField           matlab.ui.control.NumericEditField
        SamplingRateEditFieldLabel      matlab.ui.control.Label
        ResetallfiltersButton           matlab.ui.control.Button
        GroupIDEditField                matlab.ui.control.NumericEditField
        GroupIDEditFieldLabel           matlab.ui.control.Label
        SDCheckBox                      matlab.ui.control.CheckBox
        SCCheckBox                      matlab.ui.control.CheckBox
        SBCheckBox                      matlab.ui.control.CheckBox
        SACheckBox                      matlab.ui.control.CheckBox
        SubtypeLabel                    matlab.ui.control.Label
        DualVVCheckBox                  matlab.ui.control.CheckBox
        DualVHCheckBox                  matlab.ui.control.CheckBox
        DualHVCheckBox                  matlab.ui.control.CheckBox
        DualHHCheckBox                  matlab.ui.control.CheckBox
        HHHVCheckBox                    matlab.ui.control.CheckBox
        VVHHCheckBox                    matlab.ui.control.CheckBox
        HHCheckBox                      matlab.ui.control.CheckBox
        VVCheckBox                      matlab.ui.control.CheckBox
        PolarizationLabel               matlab.ui.control.Label
        WVCheckBox                      matlab.ui.control.CheckBox
        S6CheckBox                      matlab.ui.control.CheckBox
        S5CheckBox                      matlab.ui.control.CheckBox
        S4CheckBox                      matlab.ui.control.CheckBox
        S3CheckBox                      matlab.ui.control.CheckBox
        S2CheckBox                      matlab.ui.control.CheckBox
        S1CheckBox                      matlab.ui.control.CheckBox
        EWCheckBox                      matlab.ui.control.CheckBox
        IWCheckBox                      matlab.ui.control.CheckBox
        BeamModeLabel                   matlab.ui.control.Label
        DESCENDINGCheckBox              matlab.ui.control.CheckBox
        ASCENDINGCheckBox               matlab.ui.control.CheckBox
        OrbitDirectionLabel             matlab.ui.control.Label
        XMLMetadataOCNCheckBox          matlab.ui.control.CheckBox
        XMLMetadataSLCCheckBox          matlab.ui.control.CheckBox
        XMLMetadataGRDMDCheckBox        matlab.ui.control.CheckBox
        XMLMetadataGRDHSCheckBox        matlab.ui.control.CheckBox
        XMLMetadataRAWCheckBox          matlab.ui.control.CheckBox
        XMLMetadataGRDHDCheckBox        matlab.ui.control.CheckBox
        XMLMetadataGRDMSCheckBox        matlab.ui.control.CheckBox
        L0RawDataRAWCheckBox            matlab.ui.control.CheckBox
        L2OceanOCNCheckBox              matlab.ui.control.CheckBox
        L1SingleLookComplexSLCCheckBox  matlab.ui.control.CheckBox
        L1DetectedHighResSinglePolGRDHSCheckBox  matlab.ui.control.CheckBox
        L1DetectedMidResSinglePolGRDMSCheckBox  matlab.ui.control.CheckBox
        L1DetectedMidResDualPolGRDMDCheckBox  matlab.ui.control.CheckBox
        L1DetectedHighResDualPolGRDHDCheckBox  matlab.ui.control.CheckBox
        FileTypeLabel                   matlab.ui.control.Label
        ClearPathFrameButton            matlab.ui.control.Button
        EnddateDatePicker               matlab.ui.control.DatePicker
        EnddateDatePickerLabel          matlab.ui.control.Label
        StartdateDatePicker             matlab.ui.control.DatePicker
        StartdateDatePickerLabel        matlab.ui.control.Label
        DatasetDropDown                 matlab.ui.control.DropDown
        DatasetDropDownLabel            matlab.ui.control.Label
        FrameEndEditField               matlab.ui.control.NumericEditField
        FrameEndEditFieldLabel          matlab.ui.control.Label
        FrameStartEditField             matlab.ui.control.NumericEditField
        FrameStartEditFieldLabel        matlab.ui.control.Label
        PathEndEditField                matlab.ui.control.NumericEditField
        PathEndEditFieldLabel           matlab.ui.control.Label
        PathStartEditField              matlab.ui.control.NumericEditField
        PathStartEditFieldLabel         matlab.ui.control.Label
        DownloaderMapPanel              matlab.ui.container.Panel
        GlobalVariablesTab_3            matlab.ui.container.Tab
        PythonEnvironmentDropDown       matlab.ui.control.DropDown
        PythonEnvironmentDropDownLabel  matlab.ui.control.Label
        Label_2                         matlab.ui.control.Label
        CustomPythonEnvironmentEditField  matlab.ui.control.EditField
        UpdateAlreadyProcessedDataCheckBox  matlab.ui.control.CheckBox
        UpdateSearchFromLabel           matlab.ui.control.Label
        UpdateSearchEndDateDatePicker   matlab.ui.control.DatePicker
        SearchUpdateImagesButton        matlab.ui.control.Button
        DownloadUpdateImagesButton      matlab.ui.control.Button
        DownloadRunUpdateImagesButton   matlab.ui.control.Button
        UpdateSearchStatusLabel         matlab.ui.control.Label
        UpdateSearchParametersLabel     matlab.ui.control.Label
        UpdateSearchParametersTextArea  matlab.ui.control.TextArea
        UpdateImagesTableLabel          matlab.ui.control.Label
        SelectAllUpdateImagesButton     matlab.ui.control.Button
        DeselectAllUpdateImagesButton   matlab.ui.control.Button
        UpdateImagesTable               matlab.ui.control.Table
        PythonEnvironmentLabel          matlab.ui.control.Label
        AOITab_2                        matlab.ui.container.Tab
        SENMapPanel                     matlab.ui.container.Panel
        ResetSENMapButton               matlab.ui.control.Button
        AOIboxboundariesLabel_3         matlab.ui.control.Label
        MaxlatitudeEditField            matlab.ui.control.NumericEditField
        MaxlatitudeEditFieldLabel       matlab.ui.control.Label
        MinlatitudeEditField            matlab.ui.control.NumericEditField
        MinlatitudeEditFieldLabel       matlab.ui.control.Label
        MaxlongitudeEditField           matlab.ui.control.NumericEditField
        MaxlongitudeEditFieldLabel      matlab.ui.control.Label
        MinlongitudeEditField           matlab.ui.control.NumericEditField
        MinlongitudeEditFieldLabel      matlab.ui.control.Label
        MasterProcessingTab             matlab.ui.container.Tab
        AutoMasterCheckBox              matlab.ui.control.CheckBox
        PolarisationDropDown            matlab.ui.control.DropDown
        PolarisationDropDownLabel       matlab.ui.control.Label
        MasterprocessingCheckBox        matlab.ui.control.CheckBox
        MasterdateDatePicker            matlab.ui.control.DatePicker
        MasterdateDatePickerLabel       matlab.ui.control.Label
        SlavesProcessingTab             matlab.ui.container.Tab
        BrowseButton_2                  matlab.ui.control.Button
        BrowseButton                    matlab.ui.control.Button
        FirststepDropDown               matlab.ui.control.DropDown
        FirststepDropDownLabel          matlab.ui.control.Label
        onlyforExternalDEMLabel_2       matlab.ui.control.Label
        onlyforExternalDEMLabel         matlab.ui.control.Label
        DEMresamplingmethodDropDown     matlab.ui.control.DropDown
        DEMresamplingmethodDropDownLabel  matlab.ui.control.Label
        DEMcoregpathEditField           matlab.ui.control.EditField
        DEMcoregpathEditFieldLabel      matlab.ui.control.Label
        DEMifgpathEditField             matlab.ui.control.EditField
        DEMifgpathEditFieldLabel        matlab.ui.control.Label
        DEMcoregistrationDropDown       matlab.ui.control.DropDown
        DEMcoregistrationDropDownLabel  matlab.ui.control.Label
        DEMinterferogramDropDown        matlab.ui.control.DropDown
        DEMinterferogramDropDownLabel   matlab.ui.control.Label
        SlavesremovalafterprocessingCheckBox  matlab.ui.control.CheckBox
        CoherenceandLIATab              matlab.ui.container.Tab
        Label_3                         matlab.ui.control.Label
        EPSGcodeEditField               matlab.ui.control.NumericEditField
        EPSGcodeEditFieldLabel          matlab.ui.control.Label
        TerraincorrectedCoherenceandLIACheckBox  matlab.ui.control.CheckBox
        ComputationalResourcesTab       matlab.ui.container.Tab
        CacheEditField                  matlab.ui.control.EditField
        CacheEditFieldLabel             matlab.ui.control.Label
        RAMtobeusedintheprocessingGBformatnnGLabel  matlab.ui.control.Label
        CPUEditField                    matlab.ui.control.NumericEditField
        CPUEditFieldLabel               matlab.ui.control.Label
        NumberofcorestobeusedintheprocessingLabel  matlab.ui.control.Label
        PathEditField                   matlab.ui.control.EditField
        PathEditFieldLabel              matlab.ui.control.Label
        FullpathtoSNAPgptfolderLabel    matlab.ui.control.Label
        ImagesTab_SEN                   matlab.ui.container.Tab
        OpenSENSlavesFolderButton       matlab.ui.control.Button
        ImportSENButton                 matlab.ui.control.Button
        Label_SENImages                 matlab.ui.control.Label
        ImportedSENImagesTable          matlab.ui.control.Table
        SaveLoadTab                     matlab.ui.container.Tab
        StatusLamp_3                    matlab.ui.control.Lamp
        StatusLamp_3Label               matlab.ui.control.Label
        LoadButton                      matlab.ui.control.Button
        SavetheconfiguredparametersfortheInSARpreprocessingLabel_3  matlab.ui.control.Label
        StatusLamp                      matlab.ui.control.Lamp
        StatusLampLabel                 matlab.ui.control.Label
        SavetheconfiguredparametersfortheInSARpreprocessingLabel  matlab.ui.control.Label
        SaveButton                      matlab.ui.control.Button
        RunTab                          matlab.ui.container.Tab
        ThecodewillstopattheendofthecurrentstepLabel_2  matlab.ui.control.Label
        StopButton                      matlab.ui.control.Button
        MessagesTextArea                matlab.ui.control.TextArea
        MessagesTextArea_3Label         matlab.ui.control.Label
        PreprocessingstatusLamp         matlab.ui.control.Lamp
        PreprocessingstatusLampLabel    matlab.ui.control.Label
        Label_7                         matlab.ui.control.Label
        StartButton                     matlab.ui.control.Button
        Label_6                         matlab.ui.control.Label
        CosmoSkyMedPanel                matlab.ui.container.Panel
        TabGroup2                       matlab.ui.container.TabGroup
        GlobalVariablesTab_4            matlab.ui.container.Tab
        PythonEnvironmentDropDown_2     matlab.ui.control.DropDown
        PythonEnvironmentDropDown_2Label  matlab.ui.control.Label
        Label_4                         matlab.ui.control.Label
        CustomPythonEnvironmentEditField_2  matlab.ui.control.EditField
        CustomPythonEnvironmentLabel    matlab.ui.control.Label
        AOITab                          matlab.ui.container.Tab
        AOIboxboundariesLabel_2         matlab.ui.control.Label
        MinlatitudeEditField_2          matlab.ui.control.NumericEditField
        MinlatitudeEditField_2Label     matlab.ui.control.Label
        MaxlatitudeEditField_2          matlab.ui.control.NumericEditField
        MaxlatitudeEditField_2Label     matlab.ui.control.Label
        MaxlongitudeEditField_2         matlab.ui.control.NumericEditField
        MaxlongitudeEditField_2Label    matlab.ui.control.Label
        MinlongitudeEditField_2         matlab.ui.control.NumericEditField
        MinlongitudeEditField_2Label    matlab.ui.control.Label
        DrawCSKAOIButton                matlab.ui.control.Button
        Panel                           matlab.ui.container.Panel
        MasterProcessingTab_2           matlab.ui.container.Tab
        AutoMasterCheckBox_2            matlab.ui.control.CheckBox
        MasterprocessingCheckBox_2      matlab.ui.control.CheckBox
        MasterdateDatePicker_2          matlab.ui.control.DatePicker
        MasterdateDatePicker_2Label     matlab.ui.control.Label
        SlavesProcessingTab_2           matlab.ui.container.Tab
        BrowseButton_3                  matlab.ui.control.Button
        CoregistrationGCPsnumberEditField  matlab.ui.control.NumericEditField
        CoregistrationGCPsnumberEditFieldLabel  matlab.ui.control.Label
        FirststepDropDown_2             matlab.ui.control.DropDown
        FirststepDropDown_2Label        matlab.ui.control.Label
        onlyforExternalDEMLabel_3       matlab.ui.control.Label
        DEMifgpathEditField_2           matlab.ui.control.EditField
        DEMifgpathEditField_2Label      matlab.ui.control.Label
        DEMinterferogramDropDown_2      matlab.ui.control.DropDown
        DEMinterferogramDropDown_2Label  matlab.ui.control.Label
        SlavesremovalafterprocessingCheckBox_2  matlab.ui.control.CheckBox
        CoherenceandLIATab_2            matlab.ui.container.Tab
        Label_5                         matlab.ui.control.Label
        EPSGcodeEditField_2             matlab.ui.control.NumericEditField
        EPSGcodeEditField_2Label        matlab.ui.control.Label
        TerraincorrectedCoherenceandLIACheckBox_2  matlab.ui.control.CheckBox
        ComputationalResourcesTab_2     matlab.ui.container.Tab
        CacheEditField_2                matlab.ui.control.EditField
        CacheEditField_2Label           matlab.ui.control.Label
        RAMtobeusedintheprocessingGBformatnnGLabel_2  matlab.ui.control.Label
        CPUEditField_2                  matlab.ui.control.NumericEditField
        CPUEditField_2Label             matlab.ui.control.Label
        NumberofcorestobeusedintheprocessingLabel_2  matlab.ui.control.Label
        PathEditField_2                 matlab.ui.control.EditField
        PathEditField_2Label            matlab.ui.control.Label
        FullpathtoSNAPgptfolderLabel_2  matlab.ui.control.Label
        ImagesTab_CSK                   matlab.ui.container.Tab
        OpenslavesfolderButton          matlab.ui.control.Button
        ImportCSKButton                 matlab.ui.control.Button
        Label_10                        matlab.ui.control.Label
        ImportedImagesTable             matlab.ui.control.Table
        SaveLoadTab_2                   matlab.ui.container.Tab
        StatusLamp_4                    matlab.ui.control.Lamp
        StatusLamp_4Label               matlab.ui.control.Label
        SavetheconfiguredparametersfortheInSARpreprocessingLabel_4  matlab.ui.control.Label
        LoadButton_2                    matlab.ui.control.Button
        StatusLamp_2                    matlab.ui.control.Lamp
        StatusLamp_2Label               matlab.ui.control.Label
        SavetheconfiguredparametersfortheInSARpreprocessingLabel_2  matlab.ui.control.Label
        SaveButton_2                    matlab.ui.control.Button
        RunTab_2                        matlab.ui.container.Tab
        ThecodewillstopattheendofthecurrentstepLabel  matlab.ui.control.Label
        StopButton_2                    matlab.ui.control.Button
        MessagesTextArea_2              matlab.ui.control.TextArea
        MessagesTextAreaLabel           matlab.ui.control.Label
        PreprocessingstatusLamp_2       matlab.ui.control.Lamp
        PreprocessingstatusLamp_2Label  matlab.ui.control.Label
        Label_9                         matlab.ui.control.Label
        StartButton_2                   matlab.ui.control.Button
        Label_8                         matlab.ui.control.Label
        ConstellationSwitchLabel        matlab.ui.control.Label
        ContextMenu                     matlab.ui.container.ContextMenu
        Menu                            matlab.ui.container.Menu
        Menu2                           matlab.ui.container.Menu
        ExternalLogCallback              = []
        ExternalProgressCallback         = []
    end


    properties (Access = public)

        % Variables initialization with a default value
        constellation = 'SEN';
        StopFlag = false;

        roi_CSK  % Stores the red box for Cosmo-SkyMed
        roi_SEN  % Stores the red box for Sentinel-1

        % Sentinel-1 case
        os_SEN = 0;
        python_SEN = 'python';
        update_processed_data_SEN = 0;
        update_latest_date_SEN = '';
        update_reference_zip_SEN = '';
        update_local_dates_SEN = {};
        update_search_results_SEN = struct([]);
        master_date_SEN = '20200722';
        auto_master_SEN = 1;
        master_processing_SEN = 0;
        polarisation_SEN = 'VV';
        lon_min_SEN = -180.000;
        lat_min_SEN = -90.000;
        lon_max_SEN = 180.000;
        lat_max_SEN = 90.000;
        slaves_removal_SEN = 1;
        dem_name_SEN = 'SRTM 1Sec HGT';
        dem_file_SEN = '';
        dem_name_coreg_SEN = 'SRTM 1Sec HGT';
        dem_file_coreg_SEN = '';
        dem_resampling_SEN = 'NEAREST_NEIGHBOUR';
        first_step_SEN = 1;
        coherence_tc_SEN = 0;
        epsg_code_SEN = 32633;
        gptbin_path_SEN = 'C:\Program Files\snap\bin\gpt';
        cpu_SEN = 8;
        cache_SEN = '26G';

        % Cosmo-SkyMed case
        os_CSK = 0;
        python_CSK = 'python';
        master_date_CSK = '20200722';
        auto_master_CSK = 1;
        master_processing_CSK = 0;
        lon_min_CSK = -180.000;
        lat_min_CSK = -90.000;
        lon_max_CSK = 180.000;
        lat_max_CSK = 90.000;
        slaves_removal_CSK = 1;
        dem_name_CSK = 'SRTM 1Sec HGT';
        dem_file_CSK = '';
        num_gcp_CSK = 10000;
        first_step_CSK = 1;
        coherence_tc_CSK = 0;
        epsg_code_CSK = 32633;
        gptbin_path_CSK = 'C:\Program Files\snap\bin\gpt';
        cpu_CSK = 8;
        cache_CSK = '26G';

        % Downloader tab
        roi_Downloader
        DownloaderAOIType = ""
        DownloaderAOICorners = []
        DownloaderPolygonCoords = []
        DownloaderGeoAxes
        SearchPreviewProducts = []
        PreviewFootprintShapes = []
        SelectedFootprintShapes = []
        FootprintsVisible = false

    end

    methods (Access = public)

        function updateOutput(app, message)   % Function to mirror the command window

            if app.constellation == "SEN"

                % Append the message to the existing text in the TextArea
                currentText = app.MessagesTextArea.Value;
                newText = [currentText; message];
                app.MessagesTextArea.Value = newText;

            elseif app.constellation == "CSK"

                % Append the message to the existing text in the TextArea
                currentText = app.MessagesTextArea_2.Value;
                newText = [currentText; message];
                app.MessagesTextArea_2.Value = newText;

            end

            if ~isempty(app.ExternalLogCallback)
                try
                    app.ExternalLogCallback(message);
                catch callbackError
                    warning('PHASE:PreprocessingBetaLogCallback', ...
                        'Could not forward engine log: %s', callbackError.message);
                end
            end

        end

        % Downloader parameters

        function params = getDownloaderParams(app)

            params.repository = "ASF";
            params.dataset = app.DatasetDropDown.Value;

            params.processingLevel = getSelectedFileTypes(app);
            params.beamMode = getSelectedBeamModes(app);
            params.polarization = getSelectedPolarizations(app);
            params.flightDirection = getSelectedOrbitDirections(app);
            params.subtype = getSelectedSubtypes(app);

            params.startDate = optionalDate(app, app.StartdateDatePicker.Value);
            params.endDate = optionalDate(app, app.EnddateDatePicker.Value);

            params.pathStart = optionalNumber(app, app.PathStartEditField.Value);
            params.pathEnd = optionalNumber(app, app.PathEndEditField.Value);
            params.frameStart = optionalNumber(app, app.FrameStartEditField.Value);
            params.frameEnd = optionalNumber(app, app.FrameEndEditField.Value);
            params.groupID = optionalNumber(app, app.GroupIDEditField.Value);
            params.sampling.rate = optionalNumber(app, app.SamplingRateEditField.Value);

            if app.DownloaderAOIType == "rectangle"
                params.aoi.type = "rectangle";
                params.aoi.corners = app.DownloaderAOICorners;
            elseif app.DownloaderAOIType == "polygon"
                params.aoi.type = "polygon";
                params.aoi.corners = app.DownloaderPolygonCoords;
            else
                params.aoi.type = "";
                params.aoi.corners = [];
            end


            if isempty(params.sampling.rate)
                params.sampling.unit = [];
            else
                params.sampling.unit = string(app.SamplingUnitDropDown.Value);
            end

        end

        % Get filetype

        function selected = getSelectedFileTypes(app)

            selected = {};

            if app.L1DetectedHighResDualPolGRDHDCheckBox.Value
                selected{end+1} = "GRD-HD";
            end
            if app.L1DetectedMidResDualPolGRDMDCheckBox.Value
                selected{end+1} = "GRD-MD";
            end
            if app.L1DetectedMidResSinglePolGRDMSCheckBox.Value
                selected{end+1} = "GRD-MS";
            end
            if app.L1DetectedHighResSinglePolGRDHSCheckBox.Value
                selected{end+1} = "GRD-HS";
            end
            if app.L1SingleLookComplexSLCCheckBox.Value
                selected{end+1} = "SLC";
            end
            if app.L2OceanOCNCheckBox.Value
                selected{end+1} = "OCN";
            end
            if app.L0RawDataRAWCheckBox.Value
                selected{end+1} = "RAW";
            end
            if app.XMLMetadataGRDMSCheckBox.Value
                selected{end+1} = "XML_GRD-MS";
            end
            if app.XMLMetadataGRDHDCheckBox.Value
                selected{end+1} = "XML_GRD-HD";
            end
            if app.XMLMetadataRAWCheckBox.Value
                selected{end+1} = "XML_RAW";
            end
            if app.XMLMetadataGRDHSCheckBox.Value
                selected{end+1} = "XML_GRD-HS";
            end
            if app.XMLMetadataGRDMDCheckBox.Value
                selected{end+1} = "XML_GRD-MD";
            end
            if app.XMLMetadataSLCCheckBox.Value
                selected{end+1} = "XML_SLC";
            end
            if app.XMLMetadataOCNCheckBox.Value
                selected{end+1} = "XML_OCN";
            end

        end

        % Get beam mode

        function selected = getSelectedBeamModes(app)

            selected = {};

            if app.IWCheckBox.Value
                selected{end+1} = "IW";
            end
            if app.EWCheckBox.Value
                selected{end+1} = "EW";
            end
            if app.S1CheckBox.Value
                selected{end+1} = "S1";
            end
            if app.S2CheckBox.Value
                selected{end+1} = "S2";
            end
            if app.S3CheckBox.Value
                selected{end+1} = "S3";
            end
            if app.S4CheckBox.Value
                selected{end+1} = "S4";
            end
            if app.S5CheckBox.Value
                selected{end+1} = "S5";
            end
            if app.S6CheckBox.Value
                selected{end+1} = "S6";
            end
            if app.WVCheckBox.Value
                selected{end+1} = "WV";
            end

        end

        % Get polarization

        function selected = getSelectedPolarizations(app)

            selected = {};

            if app.VVCheckBox.Value
                selected{end+1} = "VV";
            end
            if app.HHCheckBox.Value
                selected{end+1} = "HH";
            end
            if app.VVHHCheckBox.Value
                selected{end+1} = "VV+HH";
            end
            if app.HHHVCheckBox.Value
                selected{end+1} = "HH+HV";
            end
            if app.DualHHCheckBox.Value
                selected{end+1} = "Dual HH";
            end
            if app.DualHVCheckBox.Value
                selected{end+1} = "Dual HV";
            end
            if app.DualVHCheckBox.Value
                selected{end+1} = "Dual VH";
            end
            if app.DualVVCheckBox.Value
                selected{end+1} = "Dual VV";
            end

        end

        % Get orbit direction

        function selected = getSelectedOrbitDirections(app)

            selected = {};

            if app.ASCENDINGCheckBox.Value
                selected{end+1} = "ASCENDING";
            end
            if app.DESCENDINGCheckBox.Value
                selected{end+1} = "DESCENDING";
            end

        end

        %Get subtype

        function selected = getSelectedSubtypes(app)

            selected = {};

            if app.SACheckBox.Value
                selected{end+1} = "SA";
            end
            if app.SBCheckBox.Value
                selected{end+1} = "SB";
            end
            if app.SCCheckBox.Value
                selected{end+1} = "SC";
            end
            if app.SDCheckBox.Value
                selected{end+1} = "SD";
            end

        end

        function value = optionalNumber(app, inputValue)

            if isempty(inputValue) || inputValue == 0
                value = [];
            else
                value = inputValue;
            end

        end

        function value = optionalDate(app, inputValue)

            if isempty(inputValue) || ismissing(inputValue)
                value = [];
            else
                value = string(inputValue, "yyyy-MM-dd");
            end

        end

        % clear all filters
        function resetAllFilters(app)

            app.DatasetDropDown.Value = "SENTINEL-1";

            % file type
            app.L1DetectedHighResDualPolGRDHDCheckBox.Value = false;
            app.L1DetectedMidResDualPolGRDMDCheckBox.Value = false;
            app.L1DetectedMidResSinglePolGRDMSCheckBox.Value = false;
            app.L1DetectedHighResSinglePolGRDHSCheckBox.Value = false;
            app.L1SingleLookComplexSLCCheckBox.Value = false;
            app.L2OceanOCNCheckBox.Value = false;
            app.L0RawDataRAWCheckBox.Value = false;

            app.XMLMetadataGRDMSCheckBox.Value = false;
            app.XMLMetadataGRDHDCheckBox.Value = false;
            app.XMLMetadataRAWCheckBox.Value = false;
            app.XMLMetadataGRDHSCheckBox.Value = false;
            app.XMLMetadataGRDMDCheckBox.Value = false;
            app.XMLMetadataSLCCheckBox.Value = false;
            app.XMLMetadataOCNCheckBox.Value = false;

            % beam mode
            app.IWCheckBox.Value = false;
            app.EWCheckBox.Value = false;
            app.S1CheckBox.Value = false;
            app.S2CheckBox.Value = false;
            app.S3CheckBox.Value = false;
            app.S4CheckBox.Value = false;
            app.S5CheckBox.Value = false;
            app.S6CheckBox.Value = false;
            app.WVCheckBox.Value = false;

            % polarization
            app.VVCheckBox.Value = false;
            app.HHCheckBox.Value = false;
            app.VVHHCheckBox.Value = false;
            app.HHHVCheckBox.Value = false;
            app.DualHHCheckBox.Value = false;
            app.DualHVCheckBox.Value = false;
            app.DualVHCheckBox.Value = false;
            app.DualVVCheckBox.Value = false;

            % orbit
            app.ASCENDINGCheckBox.Value = false;
            app.DESCENDINGCheckBox.Value = false;

            % subtype
            app.SACheckBox.Value = false;
            app.SBCheckBox.Value = false;
            app.SCCheckBox.Value = false;
            app.SDCheckBox.Value = false;

            app.StartdateDatePicker.Value = NaT;
            app.EnddateDatePicker.Value = NaT;

            app.PathStartEditField.Value = 0;
            app.PathEndEditField.Value = 0;
            app.FrameStartEditField.Value = 0;
            app.FrameEndEditField.Value = 0;

            % sampling rate
            app.SamplingRateEditField.Value = 0;
            app.SamplingUnitDropDown.Value = "Month";

        end
        % Function to start Parameter filter with default values
        function initializeParams(app)


            app.DatasetDropDown.Value = "SENTINEL-1";

            % file type
            app.L1DetectedHighResDualPolGRDHDCheckBox.Value = false;
            app.L1DetectedMidResDualPolGRDMDCheckBox.Value = false;
            app.L1DetectedMidResSinglePolGRDMSCheckBox.Value = false;
            app.L1DetectedHighResSinglePolGRDHSCheckBox.Value = false;
            app.L1SingleLookComplexSLCCheckBox.Value = true;
            app.L2OceanOCNCheckBox.Value = false;
            app.L0RawDataRAWCheckBox.Value = false;

            app.XMLMetadataGRDMSCheckBox.Value = false;
            app.XMLMetadataGRDHDCheckBox.Value = false;
            app.XMLMetadataRAWCheckBox.Value = false;
            app.XMLMetadataGRDHSCheckBox.Value = false;
            app.XMLMetadataGRDMDCheckBox.Value = false;
            app.XMLMetadataSLCCheckBox.Value = false;
            app.XMLMetadataOCNCheckBox.Value = false;

            % beam mode
            app.IWCheckBox.Value = true;
            app.EWCheckBox.Value = false;
            app.S1CheckBox.Value = false;
            app.S2CheckBox.Value = false;
            app.S3CheckBox.Value = false;
            app.S4CheckBox.Value = false;
            app.S5CheckBox.Value = false;
            app.S6CheckBox.Value = false;
            app.WVCheckBox.Value = false;

            % polarization
            app.VVCheckBox.Value = false;
            app.HHCheckBox.Value = false;
            app.VVHHCheckBox.Value = false;
            app.HHHVCheckBox.Value = false;
            app.DualHHCheckBox.Value = false;
            app.DualHVCheckBox.Value = false;
            app.DualVHCheckBox.Value = false;
            app.DualVVCheckBox.Value = false;

            % orbit
            app.ASCENDINGCheckBox.Value = false;
            app.DESCENDINGCheckBox.Value = false;

            % subtype
            app.SACheckBox.Value = false;
            app.SBCheckBox.Value = false;
            app.SCCheckBox.Value = false;
            app.SDCheckBox.Value = false;

            app.StartdateDatePicker.Value = NaT;
            app.EnddateDatePicker.Value = NaT;

            app.PathStartEditField.Value = 0;
            app.PathEndEditField.Value = 0;
            app.FrameStartEditField.Value = 0;
            app.FrameEndEditField.Value = 0;

            % sampling rate
            app.SamplingRateEditField.Value = 0;
            app.SamplingUnitDropDown.Value = "Month";

        end

        % Function to initialize map in Downloader tab

        function initializeDownloaderMap(app)

            delete(app.DownloaderMapPanel.Children);

            app.DownloaderGeoAxes = geoaxes(app.DownloaderMapPanel);
            gx = app.DownloaderGeoAxes;

            geobasemap(gx, "satellite");
            gx.Interactions = [panInteraction; zoomInteraction];

            geolimits(gx, [35 72], [-25 45]);

        end

        function updateDownloaderRectangleCoords(app, roi)

            pos = roi.Position;

            % Position format used by PHASE map: [lat lon height width]
            minLat = pos(1);
            minLon = pos(2);
            maxLat = pos(1) + pos(3);
            maxLon = pos(2) + pos(4);

            app.MinLongitudeEditField.Value = minLon;
            app.MaxLongitudeEditField.Value = maxLon;
            app.MinLatitudeEditField.Value = minLat;
            app.MaxLatitudeEditField.Value = maxLat;

            app.DownloaderAOIType = "rectangle";

            app.DownloaderAOICorners = [
                minLon minLat;
                maxLon minLat;
                maxLon maxLat;
                minLon maxLat
            ];

            app.DownloaderPolygonCoords = [];

        end

        function updateDownloaderPolygonCoords(app, roi)

            pos = roi.Position;


            lats = pos(:,1);
            lons = pos(:,2);

            app.DownloaderAOIType = "polygon";
            app.DownloaderPolygonCoords = [lons lats];

            app.DownloaderAOICorners = [];

        end

        function AOIBoundaryValueChanged(app, event)

            minLon = app.MinLongitudeEditField.Value;
            maxLon = app.MaxLongitudeEditField.Value;
            minLat = app.MinLatitudeEditField.Value;
            maxLat = app.MaxLatitudeEditField.Value;

            if minLon >= maxLon || minLat >= maxLat
                disp("Invalid AOI boundaries.");
                return;
            end

            clearDownloaderAOIShape(app);


            app.roi_Downloader = drawrectangle(app.DownloaderGeoAxes, ...
                'Position', [minLat, minLon, maxLat-minLat, maxLon-minLon], ...
                'Color', 'r', ...
                'FaceAlpha', 0.2);

            updateDownloaderRectangleCoords(app, app.roi_Downloader);

            addlistener(app.roi_Downloader, ...
                'ROIMoved', ...
                @(src, event) updateDownloaderRectangleCoords(app, src));

        end

        function clearDownloaderAOIShape(app)

            if ~isempty(app.roi_Downloader) && isvalid(app.roi_Downloader)
                delete(app.roi_Downloader);
            end

            app.roi_Downloader = [];

        end

        %check if aoi is selected
        function hasAOI = hasDownloaderAOI(app)

            hasAOI = false;

            if app.DownloaderAOIType == "rectangle" && ~isempty(app.DownloaderAOICorners)
                hasAOI = true;
                return;
            end

            if app.DownloaderAOIType == "polygon" && ~isempty(app.DownloaderPolygonCoords)
                hasAOI = true;
                return;
            end

        end
        % Download confirmation popup
        function confirmed = showDownloadConfirmation(app)

            summaryFile = fullfile(pwd, "downloadasf", "search_summary.json");

            if ~exist(summaryFile, "file")
                uialert(app.UIFigure, ...
                    "Could not find search_summary.json.", ...
                    "Search Summary Missing", ...
                    "Icon", "error");
                confirmed = false;
                return;
            end

            summary = jsondecode(fileread(summaryFile));

            imageCount = summary.product_count;
            totalSizeGB = summary.total_size_gb;

            message = sprintf( ...
            "ASF found %d SAR images.%s%sEstimated total download size: %.2f GB%s%sDo you want to continue with the download?", ...
            imageCount, newline, newline, totalSizeGB, newline, newline);

            answer = uiconfirm(app.UIFigure, ...
                message, ...
                "Confirm Download", ...
                "Options", ["Download", "Cancel"], ...
                "DefaultOption", "Download", ...
                "CancelOption", "Cancel");

            confirmed = strcmp(answer, "Download");

        end

        % Preview confirmation popup
        function previewRequested = showSearchPreviewChoice(app)

            summaryFile = fullfile(pwd, "downloadasf", "search_summary.json");

            if ~exist(summaryFile, "file")
                uialert(app.UIFigure, ...
                    "Could not find search_summary.json.", ...
                    "Search Summary Missing", ...
                    "Icon", "error");

                previewRequested = false;
                return;
            end

            summary = jsondecode(fileread(summaryFile));

            imageCount = summary.product_count;
            totalSizeGB = summary.total_size_gb;

            message = sprintf( ...
                "ASF found %d SAR images.%s%sEstimated total size: %.2f GB", ...
                imageCount, newline, newline, totalSizeGB);

            answer = uiconfirm(app.UIFigure, ...
                message, ...
                "Search Results", ...
                "Options", ["Preview Images", "OK"], ...
                "DefaultOption", "Preview Images", ...
                "CancelOption", "OK");

            previewRequested = strcmp(answer, "Preview Images");

        end

        % save last download parameters to use later
        function saveLastDownloadParams(app)

            params = getDownloaderParams(app);

            summaryFile = fullfile(pwd, "downloadasf", "search_summary.json");

            if exist(summaryFile, "file")

                summary = jsondecode(fileread(summaryFile));

                if isfield(summary, "products") && ~isempty(summary.products)

                    products = summary.products;

                    paths = unique([products.pathNumber]);
                    frames = unique([products.frameNumber]);
                    directions = unique(string({products.flightDirection}));

                    if numel(paths) == 1
                        params.pathStart = paths(1);
                        params.pathEnd = paths(1);
                    end

                    if numel(frames) == 1
                        params.frameStart = frames(1);
                        params.frameEnd = frames(1);
                    end

                    if numel(directions) == 1
                        params.flightDirection = directions;
                    end

                end

            end

            outputFolder = fullfile(pwd, "downloadasf");

            if ~exist(outputFolder, "dir")
                mkdir(outputFolder);
            end

            outputFile = fullfile(outputFolder, "last_download_request.json");

            fid = fopen(outputFile, "w");

            if fid == -1
                disp("Could not save last download parameters.");
                return;
            end

            fprintf(fid, "%s", jsonencode(params, PrettyPrint=true));
            fclose(fid);

            disp("Last download parameters saved.");

        end

        % updating the parameters to match last download's
        function applyDownloaderParams(app, params)

            resetAllFilters(app);

            app.DatasetDropDown.Value = string(params.dataset);

            setFileTypeCheckboxes(app, params.processingLevel);
            setBeamModeCheckboxes(app, params.beamMode);
            setPolarizationCheckboxes(app, params.polarization);
            setOrbitDirectionCheckboxes(app, params.flightDirection);
            setSubtypeCheckboxes(app, params.subtype);


            if ~isempty(params.startDate)
                app.StartdateDatePicker.Value = datetime(params.startDate);
            end

            if ~isempty(params.endDate)
                app.EnddateDatePicker.Value = datetime(params.endDate);
            end

            if ~isempty(params.pathStart)
                app.PathStartEditField.Value = params.pathStart;
            end

            if ~isempty(params.pathEnd)
                app.PathEndEditField.Value = params.pathEnd;
            end

            if ~isempty(params.frameStart)
                app.FrameStartEditField.Value = params.frameStart;
            end

            if ~isempty(params.frameEnd)
                app.FrameEndEditField.Value = params.frameEnd;
            end
            if isfield(params, "sampling") && ~isempty(params.sampling.rate)
                app.SamplingRateEditField.Value = params.sampling.rate;
                app.SamplingUnitDropDown.Value = string(params.sampling.unit);
            end

            if isfield(params, "aoi") && isfield(params.aoi, "type")

                clearDownloaderAOIShape(app);

                if string(params.aoi.type) == "rectangle"

                    coords = params.aoi.corners;

                    minLon = min(coords(:,1));
                    maxLon = max(coords(:,1));
                    minLat = min(coords(:,2));
                    maxLat = max(coords(:,2));

                    app.roi_Downloader = drawrectangle(app.DownloaderGeoAxes, ...
                        "Position", [minLat, minLon, maxLat-minLat, maxLon-minLon], ...
                        "Color", "r", ...
                        "FaceAlpha", 0.2);

                    updateDownloaderRectangleCoords(app, app.roi_Downloader);

                    addlistener(app.roi_Downloader, ...
                        "ROIMoved", ...
                        @(src, event) updateDownloaderRectangleCoords(app, src));

                elseif string(params.aoi.type) == "polygon"

                    coords = params.aoi.coordinates;

                    app.roi_Downloader = drawpolygon(app.DownloaderGeoAxes, ...
                        "Position", [coords(:,2), coords(:,1)], ...
                        "Color", "r", ...
                        "FaceAlpha", 0.2);

                    updateDownloaderPolygonCoords(app, app.roi_Downloader);

                    addlistener(app.roi_Downloader, ...
                        "ROIMoved", ...
                        @(src, event) updateDownloaderPolygonCoords(app, src));

                end
            end

        end

        function setFileTypeCheckboxes(app, values)

            values = string(values);

            app.L1DetectedHighResDualPolGRDHDCheckBox.Value = any(values == "GRD-HD");
            app.L1DetectedMidResDualPolGRDMDCheckBox.Value = any(values == "GRD-MD");
            app.L1DetectedMidResSinglePolGRDMSCheckBox.Value = any(values == "GRD-MS");
            app.L1DetectedHighResSinglePolGRDHSCheckBox.Value = any(values == "GRD-HS");

            app.L1SingleLookComplexSLCCheckBox.Value = any(values == "SLC");
            app.L2OceanOCNCheckBox.Value = any(values == "OCN");
            app.L0RawDataRAWCheckBox.Value = any(values == "RAW");

            app.XMLMetadataGRDMSCheckBox.Value = any(values == "XML_GRD-MS");
            app.XMLMetadataGRDHDCheckBox.Value = any(values == "XML_GRD-HD");
            app.XMLMetadataRAWCheckBox.Value = any(values == "XML_RAW");
            app.XMLMetadataGRDHSCheckBox.Value = any(values == "XML_GRD-HS");
            app.XMLMetadataGRDMDCheckBox.Value = any(values == "XML_GRD-MD");
            app.XMLMetadataSLCCheckBox.Value = any(values == "XML_SLC");
            app.XMLMetadataOCNCheckBox.Value = any(values == "XML_OCN");

        end

        function setBeamModeCheckboxes(app, values)

            values = string(values);

            app.IWCheckBox.Value = any(values == "IW");
            app.EWCheckBox.Value = any(values == "EW");
            app.S1CheckBox.Value = any(values == "S1");
            app.S2CheckBox.Value = any(values == "S2");
            app.S3CheckBox.Value = any(values == "S3");
            app.S4CheckBox.Value = any(values == "S4");
            app.S5CheckBox.Value = any(values == "S5");
            app.S6CheckBox.Value = any(values == "S6");
            app.WVCheckBox.Value = any(values == "WV");

        end

        function setPolarizationCheckboxes(app, values)

            values = string(values);

            app.VVCheckBox.Value = any(values == "VV");
            app.HHCheckBox.Value = any(values == "HH");
            app.VVHHCheckBox.Value = any(values == "VV+HH");
            app.HHHVCheckBox.Value = any(values == "HH+HV");
            app.DualHHCheckBox.Value = any(values == "Dual HH");
            app.DualHVCheckBox.Value = any(values == "Dual HV");
            app.DualVHCheckBox.Value = any(values == "Dual VH");
            app.DualVVCheckBox.Value = any(values == "Dual VV");

        end

        function setOrbitDirectionCheckboxes(app, values)

            values = string(values);

            app.ASCENDINGCheckBox.Value = any(values == "ASCENDING");
            app.DESCENDINGCheckBox.Value = any(values == "DESCENDING");

        end

        function setSubtypeCheckboxes(app, values)

            values = string(values);

            app.SACheckBox.Value = any(values == "SA");
            app.SBCheckBox.Value = any(values == "SB");
            app.SCCheckBox.Value = any(values == "SC");
            app.SCCheckBox.Value = any(values == "SD");


        end

        % helper function for select all check box
        function updatePreviewSelectionState(app)


            tableData = app.PreviewResultsTable.Data;

            if isempty(tableData)

                app.SelectAllCheckBox.Value = false;
                app.DownloadSelectedButton.Enable = "off";
                app.PreviewSelectedLabel.Text = "Selected: 0";

                return;

            end

            selected = cell2mat(tableData(:,1));

            selectedCount = sum(selected);

            app.PreviewSelectedLabel.Text = ...
                "Selected: " + selectedCount;

            app.DownloadSelectedButton.Enable = ...
                matlab.lang.OnOffSwitchState(selectedCount > 0);

            app.SelectAllCheckBox.Value = all(selected);
        end

        % load search result to preview tab
        function loadSearchPreview(app)

            summaryFile = fullfile(pwd, "downloadasf", "search_summary.json");

            if ~exist(summaryFile, "file")
                uialert(app.UIFigure, ...
                    "Could not find search_summary.json.", ...
                    "Missing Results", ...
                    "Icon", "error");
                return;
            end

            summary = jsondecode(fileread(summaryFile));

            if ~isfield(summary, "products") || isempty(summary.products)
                app.PreviewResultsTable.Data = {};
                app.PreviewSelectedLabel.Text = "Selected: 0";
                app.DownloadSelectedButton.Enable = "off";
                return;
            end

            app.SearchPreviewProducts = summary.products;

            products = summary.products;
            n = numel(products);

            % show recommended path frame and direction
            paths = unique([products.pathNumber]);
            frames = unique([products.frameNumber]);
            directions = unique(string({products.flightDirection}));

            hasOnlyOneCombination = ...
                numel(paths) == 1 && ...
                numel(frames) == 1 && ...
                numel(directions) == 1;

            if hasOnlyOneCombination

                app.RecommendedFramePathLabel.Text = "";

            elseif isfield(summary, "best_path") && ...
                   isfield(summary, "best_frame") && ...
                   isfield(summary, "best_direction")

                app.RecommendedFramePathLabel.Text = sprintf( ...
                    "Recommended: Path %d, Frame %d, %s", ...
                    summary.best_path, ...
                    summary.best_frame, ...
                    char(summary.best_direction));

            else

                app.RecommendedFramePathLabel.Text = ...
                    "Recommended: No best path/frame found.";

            end

            tableData = cell(n, 7);

            for i = 1:n

                p = products(i);

                sceneName = char(p.sceneName);

                tokens = regexp(sceneName, '\d{8}T\d{6}', 'match');

                if ~isempty(tokens)
                    rawDateTime = tokens{1};

                    rawDate = rawDateTime(1:8);
                    rawTime = rawDateTime(10:15);

                    dateText = datestr(datetime(rawDate, ...
                        "InputFormat", "yyyyMMdd"), "yyyy-mm-dd");

                    timeText = [rawTime(1:2) ':' rawTime(3:4) ':' rawTime(5:6)];
                else
                    dateText = char(p.startTime);
                    timeText = '';
                end

                tableData{i,1} = false;
                tableData{i,2} = dateText;
                tableData{i,3} = timeText;
                tableData{i,4} = int32(p.pathNumber);
                tableData{i,5} = int32(p.frameNumber);
                tableData{i,6} = char(p.flightDirection);
                tableData{i,7} = char(sprintf("%.2f", double(p.size)));

            end

            app.PreviewResultsTable.ColumnName = { ...
                'Select', ...
                'Date', ...
                'Time', ...
                'Path', ...
                'Frame', ...
                'Direction', ...
                'Size GB'};

            app.PreviewResultsTable.ColumnEditable = ...
                [true false false false false false false];

            app.PreviewResultsTable.ColumnSortable = ...
                [false true true true true true true];

            app.PreviewResultsTable.Data = tableData;


        end

        function selectedProducts = getSelectedPreviewProducts(app)

            tableData = app.PreviewResultsTable.Data;

            if isempty(tableData)
                selectedProducts = [];
                return;
            end

            selected = cell2mat(tableData(:,1));
            selectedProducts = app.SearchPreviewProducts(selected);

        end

        % count existing files
        function existingCount = countExistingDownloadFiles(app, products)

            appPath = phase_preprocessing_beta.projectRoot();
            downloadFolder = fullfile(appPath, "PHASE_Preprocessing", "slaves");

            existingCount = 0;

            for i = 1:numel(products)
                filePath = fullfile(downloadFolder, string(products(i).sceneName) + ".zip");

                if exist(filePath, "file")
                    existingCount = existingCount + 1;
                end
            end

        end

        % count existing bytes in folder
        function existingBytes = countExistingDownloadBytes(app, products)

            appPath = phase_preprocessing_beta.projectRoot();
            downloadFolder = fullfile(appPath, "PHASE_Preprocessing", "slaves");

            existingBytes = 0;

            for i = 1:numel(products)
                filePath = fullfile(downloadFolder, string(products(i).sceneName) + ".zip");

                if exist(filePath, "file")
                    info = dir(filePath);
                    existingBytes = existingBytes + info.bytes;
                end
            end

        end
        % valid download
        function valid = validateDownloadCompatibility(app, products)

            valid = false;

            if isempty(products)
                uialert(app.UIFigure, ...
                    "No images selected.", ...
                    "Download Error", ...
                    "Icon", "warning");
                return;
            end

            paths = unique([products.pathNumber]);
            frames = unique([products.frameNumber]);
            directions = unique(string({products.flightDirection}));

            if numel(paths) > 1 || numel(frames) > 1 || numel(directions) > 1

                message = compose( ...
                    "Selected images are not compatible for PHASE.\n\n" + ...
                    "All selected images must have the same path, frame, and direction.\n\n" + ...
                    "Paths: %s\nFrames: %s\nDirections: %s", ...
                    mat2str(paths), ...
                    mat2str(frames), ...
                    strjoin(directions, ", "));

                uialert(app.UIFigure, message, ...
                    "Incompatible Selection", ...
                    "Icon", "error");
                return;
            end

            valid = true;

        end

        function confirmed = confirmDownloadSelection(app, products)

            fileCount = numel(products);
            totalSizeGB = sum([products.size]);

            existingCount = countExistingDownloadFiles(app, products);
            newCount = fileCount - existingCount;

            message = "Download " + string(fileCount) + " SAR images?" + newline + newline + ...
                "Total size: " + string(sprintf("%.2f", totalSizeGB)) + " GB" + newline + newline + ...
                "Already in slaves folder: " + string(existingCount) + newline + ...
                "New files to download: " + string(newCount) + newline + newline + ...
                "Files will be downloaded to:" + newline + ...
                "PHASE_Preprocessing/slaves" + newline + newline + ...
                "Old files in this folder that are not part of this download may be deleted.";

            answer = uiconfirm(app.UIFigure, ...
                char(message), ...
                "Confirm Download", ...
                "Options", ["Download", "Cancel"], ...
                "DefaultOption", "Download", ...
                "CancelOption", "Cancel");

            confirmed = strcmp(answer, "Download");

        end

        function ok = saveDownloadSelection(app, products)

            ok = false;

            downloadFolder = fullfile(pwd, "downloadasf");

            allDownloadFile = fullfile(downloadFolder, "download_data_all.json");
            selectedDownloadFile = fullfile(downloadFolder, "download_data.json");

            if ~exist(allDownloadFile, "file")
                uialert(app.UIFigure, ...
                    "Original download data was not found. Please run Search ASF again.", ...
                    "Download Error", ...
                    "Icon", "error");
                return;
            end

            allDownloadData = jsondecode(fileread(allDownloadFile));

            selectedSceneNames = string({products.sceneName});
            allInfo = allDownloadData.information;

            selectedInfo = struct([]);
            count = 0;

            for i = 1:numel(allInfo)
                if any(selectedSceneNames == string(allInfo(i).sceneName))
                    count = count + 1;
                    selectedInfo(count).sceneName = allInfo(i).sceneName;
                    selectedInfo(count).url = allInfo(i).url;
                end
            end

            output.information = selectedInfo;

            fid = fopen(selectedDownloadFile, "w");
            fprintf(fid, "%s", jsonencode(output, PrettyPrint=true));
            fclose(fid);

            % Also update search_summary.json so backend compatibility checks selected products only
            allSummaryFile = fullfile(downloadFolder, "search_summary_all.json");
            selectedSummaryFile = fullfile(downloadFolder, "search_summary.json");

            summary = jsondecode(fileread(allSummaryFile));
            summary.products = products;
            summary.product_count = numel(products);
            summary.total_size_gb = round(sum([products.size]), 2);
            summary.total_size_bytes = sum([products.size_bytes]);

            fid = fopen(selectedSummaryFile, "w");
            fprintf(fid, "%s", jsonencode(summary, PrettyPrint=true));
            fclose(fid);

            ok = true;

        end

        function ok = restoreAllDownloadFiles(app)

            ok = false;

            downloadFolder = fullfile(pwd, "downloadasf");

            allDownloadFile = fullfile(downloadFolder, "download_data_all.json");
            downloadFile = fullfile(downloadFolder, "download_data.json");

            allSummaryFile = fullfile(downloadFolder, "search_summary_all.json");
            summaryFile = fullfile(downloadFolder, "search_summary.json");

            if ~exist(allDownloadFile, "file") || ~exist(allSummaryFile, "file")
                uialert(app.UIFigure, ...
                    "Original search files were not found. Please run Search ASF again.", ...
                    "Download Error", ...
                    "Icon", "error");
                return;
            end

            copyfile(allDownloadFile, downloadFile);
            copyfile(allSummaryFile, summaryFile);

            ok = true;

        end

        function ok = waitForDeletionComplete(app, timeoutSeconds)

            ok = false;

            appPath = phase_preprocessing_beta.projectRoot();
            deletionFile = fullfile(appPath, "downloadasf", "deletion.json");

            startWait = tic;

            app.DownloadProgressLabel.Text = ...
                "Preparing download... deleting old files.";

            drawnow;

            while toc(startWait) < timeoutSeconds

                if exist(deletionFile, "file")

                    try

                        data = jsondecode(fileread(deletionFile));

                        % Backend writes plain string:
                        % "pending"
                        % "complete"

                        if ischar(data) || isstring(data)

                            if string(data) == "complete"
                                ok = true;
                                return;
                            end

                        end

                    catch
                        % File may still be written
                    end

                end

                pause(0.5);
                drawnow;

            end

            uialert(app.UIFigure, ...
                "Timed out while waiting for deletion to finish.", ...
                "Download Preparation Failed", ...
                "Icon", "error");

        end

        function runBackendDownload(app, totalBytes, products)
            disp("RUN BACKEND DOWNLOAD STARTED")
            app.DownloadProgressGauge.Visible = "on";
            app.DownloadProgressLabel.Visible = "on";
            app.DownloadProgressGauge.Value = 0;
            app.DownloadProgressLabel.Text = "Starting download...";

            appPath = phase_preprocessing_beta.projectRoot();
            controllerPath = fullfile(appPath, "downloadasf", "controller.py");

            pythonEnv = pyenv;
            pythonExe = string(pythonEnv.Executable);

            pythonExeChar = char(pythonExe);
            controllerPathChar = char(controllerPath);

            if ispc
                command = sprintf('cmd /c start "" /B "%s" "%s" download', ...
                    pythonExeChar, controllerPathChar);
            else
                command = sprintf('"%s" "%s" download &', ...
                    pythonExeChar, controllerPathChar);
            end

            disp("Download command:");
            disp(command);

            [status, cmdout] = system(command);

            disp("Download start status:");
            disp(status);
            disp("Download start output:");
            disp(cmdout);

            if status ~= 0
                app.DownloadProgressLabel.Text = "Could not start download.";
                uialert(app.UIFigure, cmdout, "Download Failed", "Icon", "error");
                return;
            end

            % waiting for backend to delete old files, before download
            if ~waitForDeletionComplete(app, 300)
                app.DownloadProgressLabel.Text = "Download preparation failed.";
                return;
            end

            existingBytesBefore = countExistingDownloadBytes(app, products);
            bytesToActuallyDownload = max(totalBytes - existingBytesBefore, 0);

            startTime = datetime("now");
            while app.DownloadProgressGauge.Value < 100
                updateDownloadProgress(app, totalBytes, startTime, existingBytesBefore, bytesToActuallyDownload);
                pause(1);
            end

            app.DownloadProgressLabel.Text = "Download completed.";

            uialert(app.UIFigure, ...
                "Download completed.", ...
                "Download Complete", ...
                "Icon", "success");

            app.initializeSENMap(); % Refresh the AOI map with the new footprints

        end

        % helper for downloading bar
        function updateDownloadProgress(app, totalBytes, startTime, existingBytesBefore, bytesToActuallyDownload)

            appPath = phase_preprocessing_beta.projectRoot();
            downloadFolder = fullfile(appPath, "PHASE_Preprocessing", "slaves");

            downloadDataFile = fullfile(appPath, "downloadasf", "download_data.json");
            downloadData = jsondecode(fileread(downloadDataFile));

            products = downloadData.information;

            currentBytes = 0;

            for i = 1:numel(products)
                filePath = fullfile(downloadFolder, string(products(i).sceneName) + ".zip");

                if exist(filePath, "file")
                    info = dir(filePath);
                    currentBytes = currentBytes + info.bytes;
                end
            end

            actualNewBytes = max(currentBytes - existingBytesBefore, 0);

            if totalBytes <= 0
                progress = 0;
            else
                progress = min(currentBytes / totalBytes, 1);
            end

            elapsedSeconds = seconds(datetime("now") - startTime);

            if elapsedSeconds > 0 && actualNewBytes > 0 && bytesToActuallyDownload > 0

                speedBytesPerSec = actualNewBytes / elapsedSeconds;
                remainingBytes = max(bytesToActuallyDownload - actualNewBytes, 0);
                remainingSeconds = remainingBytes / speedBytesPerSec;

                speedMBs = speedBytesPerSec / (1024^2);
                etaText = formatSeconds(app, remainingSeconds);

                app.DownloadProgressLabel.Text = sprintf( ...
                    "Downloaded %.1f%% | %.2f MB/s | ETA: %s", ...
                    progress * 100, speedMBs, etaText);

            elseif bytesToActuallyDownload == 0

                app.DownloadProgressLabel.Text = ...
                    "All selected files already exist in slaves.";

            else

                app.DownloadProgressLabel.Text = sprintf( ...
                    "Downloaded %.1f%% | waiting for new file download...", ...
                    progress * 100);

            end

            app.DownloadProgressGauge.Value = progress * 100;
            drawnow;

        end

        function backupSearchDownloadFiles(app)

            downloadFolder = fullfile(pwd, "downloadasf");

            filesToBackup = {
                "download_data.json", "download_data_all.json";
                "search_summary.json", "search_summary_all.json"
            };

            for i = 1:size(filesToBackup, 1)

                sourceFile = fullfile(downloadFolder, filesToBackup{i,1});
                backupFile = fullfile(downloadFolder, filesToBackup{i,2});

                if exist(sourceFile, "file")
                    copyfile(sourceFile, backupFile);
                end

            end

        end

        function text = formatSeconds(app, secondsValue)

            if isinf(secondsValue) || isnan(secondsValue)
                text = "unknown";
                return;
            end

            secondsValue = round(secondsValue);

            hours = floor(secondsValue / 3600);
            minutes = floor(mod(secondsValue, 3600) / 60);
            secondsLeft = mod(secondsValue, 60);

            if hours > 0
                text = sprintf("%dh %dm", hours, minutes);
            elseif minutes > 0
                text = sprintf("%dm %ds", minutes, secondsLeft);
            else
                text = sprintf("%ds", secondsLeft);
            end

        end

        % draw all unique footprints
        function drawAllPreviewFootprints(app)

            if isempty(app.SearchPreviewProducts)
                disp("No preview products loaded.");
                return;
            end

            latlim = app.DownloaderGeoAxes.LatitudeLimits;
            lonlim = app.DownloaderGeoAxes.LongitudeLimits;

            clearPreviewFootprints(app);

            hold(app.DownloaderGeoAxes, "on");

            products = app.SearchPreviewProducts;
            seen = strings(0);

            for i = 1:numel(products)

                p = products(i);

                key = string(p.pathNumber) + "_" + ...
                      string(p.frameNumber) + "_" + ...
                      string(p.flightDirection);

                if any(seen == key)
                    continue;
                end

                seen(end+1) = key;

                disp("Drawing footprint for:");
                disp(p.sceneName);
                disp("Path / Frame / Direction:");
                disp([p.pathNumber, p.frameNumber]);
                disp(p.flightDirection);

                [lons, lats] = getFootprintLonLat(app, p.footprint);

                disp("Plotting first coords [lat lon]:");
                disp([lats(1:min(5,end)), lons(1:min(5,end))]);

                h = geoplot(app.DownloaderGeoAxes, ...
                    lats, lons, ...
                    "y-", ...
                    "LineWidth", 2);

                app.PreviewFootprintShapes = [app.PreviewFootprintShapes; h];

                disp("Footprint plotted.");

                % TEMP TEST: zoom to first footprint
                if numel(app.PreviewFootprintShapes) == 1
                    geolimits(app.DownloaderGeoAxes, ...
                        [min(lats)-0.1, max(lats)+0.1], ...
                        [min(lons)-0.1, max(lons)+0.1]);
                end

            end

            % Comment this out while testing, because it may zoom back away.
            % geolimits(app.DownloaderGeoAxes, latlim, lonlim);
        end

        % highlight footprints for selected results
        function drawSelectedPreviewFootprints(app)

            clearSelectedFootprints(app);

            tableData = app.PreviewResultsTable.Data;

            if isempty(tableData) || isempty(app.SearchPreviewProducts)
                disp("No selected footprints to draw.");
                return;
            end

            hold(app.DownloaderGeoAxes, "on");

            selected = cell2mat(tableData(:,1));
            products = app.SearchPreviewProducts;

            seen = strings(0);

            for i = 1:numel(products)

                if ~selected(i)
                    continue;
                end

                p = products(i);

                key = string(p.pathNumber) + "_" + ...
                      string(p.frameNumber) + "_" + ...
                      string(p.flightDirection);

                if any(seen == key)
                    continue;
                end

                seen(end+1) = key;

                disp("Highlighting selected footprint for:");
                disp(p.sceneName);

                [lons, lats] = getFootprintLonLat(app, p.footprint);

                h = geoplot(app.DownloaderGeoAxes, ...
                    lats, lons, ...
                    "b-", ...
                    "LineWidth", 3);

                app.SelectedFootprintShapes = [app.SelectedFootprintShapes; h];

            end
        end

        % footprint coords helper
        function [lons, lats] = getFootprintLonLat(app, footprint)

            coords = footprint.coordinates;

            while iscell(coords)
                coords = coords{1};
            end

            % GeoJSON polygon decoded as 1 x N x 2
            if ndims(coords) == 3
                coords = squeeze(coords);
            end

            % If squeeze gives 2 x N, transpose to N x 2
            if size(coords, 1) == 2 && size(coords, 2) ~= 2
                coords = coords';
            end

            lons = coords(:,1);
            lats = coords(:,2);

            disp("Fixed footprint lon/lat min max:");
            disp([min(lons), max(lons), min(lats), max(lats)]);

        end
        % clear helpers
        function clearPreviewFootprints(app)

            for i = 1:numel(app.PreviewFootprintShapes)
                if isvalid(app.PreviewFootprintShapes(i))
                    delete(app.PreviewFootprintShapes(i));
                end
            end

            app.PreviewFootprintShapes = [];

        end

        function clearSelectedFootprints(app)

            for i = 1:numel(app.SelectedFootprintShapes)
                if isvalid(app.SelectedFootprintShapes(i))
                    delete(app.SelectedFootprintShapes(i));
                end
            end

            app.SelectedFootprintShapes = [];

        end

        function setSENUpdateControlsVisible(app, isVisible)
            if isVisible
                visibility = 'on';
            else
                visibility = 'off';
            end

            app.UpdateSearchFromLabel.Visible = visibility;
            app.UpdateSearchEndDateDatePicker.Visible = visibility;
            app.UpdateSearchParametersLabel.Visible = visibility;
            app.UpdateSearchParametersTextArea.Visible = visibility;
            app.UpdateImagesTableLabel.Visible = visibility;
            app.SelectAllUpdateImagesButton.Visible = visibility;
            app.DeselectAllUpdateImagesButton.Visible = visibility;
            app.UpdateImagesTable.Visible = visibility;
            app.SearchUpdateImagesButton.Visible = visibility;
            app.DownloadUpdateImagesButton.Visible = visibility;
            app.DownloadRunUpdateImagesButton.Visible = visibility;
            app.UpdateSearchStatusLabel.Visible = visibility;
        end

        function refreshSENUpdateContext(app)
            [latestDate, referenceZip, localDates] = app.getSENAvailableImageContext();

            app.update_reference_zip_SEN = referenceZip;
            app.update_local_dates_SEN = localDates;
            app.update_search_results_SEN = struct([]);
            app.UpdateImagesTable.Data = {};
            if isempty(app.UpdateSearchEndDateDatePicker.Value)
                app.UpdateSearchEndDateDatePicker.Value = datetime('today');
            end

            if isempty(latestDate)
                app.update_latest_date_SEN = '';
                app.UpdateSearchFromLabel.Text = 'Search new images from - to';
                app.UpdateSearchParametersTextArea.Value = {'No already processed Sentinel-1 images found.'; ...
                    'Expected data: master ZIP or slaves folders named YYYYMMDD.'};
                app.UpdateSearchStatusLabel.Text = 'Cannot search: no local Sentinel-1 image date is available.';
                app.UpdateSearchEndDateDatePicker.Enable = 'off';
                app.SearchUpdateImagesButton.Enable = 'off';
                app.DownloadUpdateImagesButton.Enable = 'off';
                app.DownloadRunUpdateImagesButton.Enable = 'off';
                return
            end

            app.update_latest_date_SEN = datestr(latestDate, 'yyyymmdd');
            app.UpdateSearchFromLabel.Text = sprintf('Search new images from %s to', datestr(latestDate, 'dd-mmm-yyyy'));
            app.UpdateSearchEndDateDatePicker.Enable = 'on';

            if isempty(referenceZip)
                app.UpdateSearchParametersTextArea.Value = {'File type: L1 Single Look Complex (SLC)'; ...
                    'Beam mode: IW'; ...
                    'Polarization: empty (all)'; ...
                    'Direction: unavailable (missing reference ZIP)'; ...
                    'Subtype: All (Sentinel-1A/B/C/D)'; ...
                    'Path Start/End: unavailable (missing reference ZIP)'; ...
                    'Frame Start/End: unavailable (missing reference ZIP)'};
                app.UpdateSearchStatusLabel.Text = 'Cannot search: master/slave Sentinel-1 ZIP reference not found.';
                app.SearchUpdateImagesButton.Enable = 'off';
                app.DownloadUpdateImagesButton.Enable = 'off';
                app.DownloadRunUpdateImagesButton.Enable = 'off';
                return
            end

            [~, refName, refExt] = fileparts(referenceZip);
            app.UpdateSearchParametersTextArea.Value = {'File type: L1 Single Look Complex (SLC)'; ...
                'Beam mode: IW'; ...
                'Polarization: empty (all)'; ...
                'Direction: same as reference image (available after search)'; ...
                'Subtype: All (Sentinel-1A/B/C/D)'; ...
                'Path Start/End: same as reference image (available after search)'; ...
                'Frame Start/End: same as reference image (available after search)'; ...
                ['Reference: ' refName refExt]};
            app.UpdateSearchStatusLabel.Text = 'Click Search Images to query ASF.';
            app.SearchUpdateImagesButton.Enable = 'on';
            app.DownloadUpdateImagesButton.Enable = 'off';
            app.DownloadRunUpdateImagesButton.Enable = 'off';
        end

        function [latestDate, referenceZip, localDates] = getSENAvailableImageContext(app)
            currentFolder = phase_preprocessing_beta.projectRoot();
            projectFolder = fullfile(currentFolder, 'PHASE_Preprocessing');
            localDates = {};
            referenceZip = '';

            masterZips = dir(fullfile(projectFolder, 'master', '*.zip'));
            for k = 1:numel(masterZips)
                zipPath = fullfile(masterZips(k).folder, masterZips(k).name);
                dateText = app.parseSENDateFromFilename(masterZips(k).name);
                if ~isempty(dateText)
                    localDates{end + 1} = dateText; %#ok<AGROW>
                end
                if isempty(referenceZip)
                    referenceZip = zipPath;
                end
            end

            slavesFolder = fullfile(projectFolder, 'slaves');
            slaveDirs = dir(slavesFolder);
            for k = 1:numel(slaveDirs)
                if slaveDirs(k).isdir && ~strcmp(slaveDirs(k).name, '.') && ~strcmp(slaveDirs(k).name, '..') && ...
                        ~isempty(regexp(slaveDirs(k).name, '^\d{8}$', 'once'))
                    localDates{end + 1} = slaveDirs(k).name; %#ok<AGROW>
                end
            end

            slaveZips = dir(fullfile(slavesFolder, '**', '*.zip'));
            for k = 1:numel(slaveZips)
                zipPath = fullfile(slaveZips(k).folder, slaveZips(k).name);
                dateText = app.parseSENDateFromFilename(slaveZips(k).name);
                if ~isempty(dateText)
                    localDates{end + 1} = dateText; %#ok<AGROW>
                end
                if isempty(referenceZip)
                    referenceZip = zipPath;
                end
            end

            localDates = unique(localDates);
            if isempty(localDates)
                latestDate = [];
            else
                dateValues = datetime(localDates, 'InputFormat', 'yyyyMMdd');
                latestDate = max(dateValues);
            end
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

        function reprojectGeoTiffIfGeographic(app, filePath, epsgCode, bandIndex)
            if nargin < 4
                bandIndex = [];
            end

            info = geotiffinfo(filePath);
            imageData = geotiffread(filePath);

            if ~isempty(bandIndex)
                if size(imageData, 3) < bandIndex
                    return
                end
                imageData = imageData(:, :, bandIndex);
            end

            imageData = double(imageData);
            R = info.SpatialRef;

            if isprop(R, 'LongitudeLimits') && isprop(R, 'LatitudeLimits')
                lonStep = (R.LongitudeLimits(2) - R.LongitudeLimits(1)) / (R.RasterSize(2) - 1);
                latStep = (R.LatitudeLimits(2) - R.LatitudeLimits(1)) / (R.RasterSize(1) - 1);
                lon = repmat(R.LongitudeLimits(1):lonStep:R.LongitudeLimits(2), R.RasterSize(1), 1);
                lat = repmat((fliplr(R.LatitudeLimits(1):latStep:R.LatitudeLimits(2)))', 1, R.RasterSize(2));
            elseif isprop(R, 'XWorldLimits') && isprop(R, 'YWorldLimits')
                geotiffwrite(filePath, imageData, R, "CoordRefSysCode", strcat('EPSG:', num2str(epsgCode)));
                return
            else
                error('Unsupported GeoTIFF spatial reference type: %s', class(R));
            end

            [x, y] = projfwd(projcrs(epsgCode), lat, lon);
            RProj = maprefcells([min(x(:)), max(x(:))], [min(y(:)), max(y(:))], size(imageData), 'ColumnsStartFrom', 'north');
            RProj.ProjectedCRS = projcrs(epsgCode);

            geotiffwrite(filePath, imageData, RProj, "CoordRefSysCode", strcat('EPSG:', num2str(epsgCode)));
        end

        function searchSENUpdateImages(app)
            app.refreshSENUpdateContext();

            if isempty(app.update_latest_date_SEN) || isempty(app.update_reference_zip_SEN)
                return
            end

            endDate = app.UpdateSearchEndDateDatePicker.Value;
            endDateText = datestr(endDate, 'yyyymmdd');
            if datetime(endDateText, 'InputFormat', 'yyyyMMdd') < datetime(app.update_latest_date_SEN, 'InputFormat', 'yyyyMMdd')
                uialert(app.UIFigure, 'The end date cannot be earlier than the latest local image date.', 'Invalid Date Range');
                return
            end

            app.UpdateSearchStatusLabel.Text = 'Searching ASF...';
            drawnow;

            currentFolder = phase_preprocessing_beta.projectRoot();
            scriptPath = fullfile(currentFolder, 'pythonScripts', 'search_update_sentinel1_images.py');
            if exist(scriptPath, 'file') ~= 2
                uialert(app.UIFigure, ['Missing helper script: ' scriptPath], 'Missing File');
                app.UpdateSearchStatusLabel.Text = 'Search failed: missing helper script.';
                return
            end

            outputFile = [tempname '.json'];
            cleanupOutput = onCleanup(@() app.deleteFileIfExists(outputFile));
            excludeDates = strjoin(app.update_local_dates_SEN, ',');
            command = app.buildPythonCommand(scriptPath, {'--reference', app.update_reference_zip_SEN, ...
                '--start-date', app.update_latest_date_SEN, '--end-date', endDateText, ...
                '--exclude-dates', excludeDates, '--output', outputFile});

            [status, commandOutput] = system(command);
            if status ~= 0
                app.UpdateSearchStatusLabel.Text = 'Search failed.';
                uialert(app.UIFigure, commandOutput, 'ASF Search Failed');
                return
            end

            if exist(outputFile, 'file') ~= 2
                app.UpdateSearchStatusLabel.Text = 'Search failed: no JSON output.';
                uialert(app.UIFigure, commandOutput, 'ASF Search Failed');
                return
            end

            searchData = jsondecode(fileread(outputFile));
            if isfield(searchData, 'params')
                app.UpdateSearchParametersTextArea.Value = app.formatSENUpdateParameters(searchData.params);
            end

            if isfield(searchData, 'results')
                app.update_search_results_SEN = searchData.results;
                app.populateSENUpdateImagesTable(searchData.results);
            else
                app.update_search_results_SEN = struct([]);
                app.UpdateImagesTable.Data = {};
                app.UpdateSearchStatusLabel.Text = 'No new matching images found.';
                app.DownloadUpdateImagesButton.Enable = 'off';
                app.DownloadRunUpdateImagesButton.Enable = 'off';
            end
        end

        function success = downloadSelectedSENUpdateImages(app, autoRun, showCompletionDialog)
            if nargin < 2
                autoRun = false;
            end
            if nargin < 3
                showCompletionDialog = true;
            end

            success = false;
            if autoRun
                app.TabGroup.SelectedTab = app.RunTab;
                drawnow;
            end

            tableData = app.UpdateImagesTable.Data;
            if isempty(tableData)
                uialert(app.UIFigure, 'No images are available for download. Run Search Images first.', 'No Images');
                return
            end

            success = app.downloadSelectedSENUpdateImagesImpl(autoRun, showCompletionDialog);
        end

        function success = downloadSelectedSENUpdateImagesImpl(app, autoRun, showCompletionDialog)
            if nargin < 2
                autoRun = false;
            end
            if nargin < 3
                showCompletionDialog = true;
            end

            success = false;
            tableData = app.UpdateImagesTable.Data;
            selected = false(size(tableData, 1), 1);
            for k = 1:size(tableData, 1)
                value = tableData{k, 1};
                if islogical(value)
                    selected(k) = value;
                elseif isnumeric(value)
                    selected(k) = value ~= 0;
                end
            end

            selectedRows = tableData(selected, :);
            if isempty(selectedRows)
                uialert(app.UIFigure, 'Select at least one image to download.', 'No Selection');
                return
            end

            currentFolder = phase_preprocessing_beta.projectRoot();
            scriptPath = fullfile(currentFolder, 'pythonScripts', 'downloader_update_images.py');
            if exist(scriptPath, 'file') ~= 2
                uialert(app.UIFigure, ['Missing helper script: ' scriptPath], 'Missing File');
                app.UpdateSearchStatusLabel.Text = 'Download failed: missing helper script.';
                updateOutput(app, 'Update download failed: missing helper script.');
                return
            end

            destinationFolder = fullfile(currentFolder, 'PHASE_Preprocessing', 'slaves');
            projectFolder = fullfile(currentFolder, 'PHASE_Preprocessing');
            if exist(destinationFolder, 'dir') ~= 7
                mkdir(destinationFolder);
            end

            urlsFile = fullfile(projectFolder, 'update_selected_urls.txt');
            outputFile = fullfile(projectFolder, 'update_download_summary.json');
            downloadProof = fullfile(destinationFolder, 'dummyUpdateDownload.txt');
            if exist(downloadProof, 'file') == 2
                delete(downloadProof);
            end
            cleanupDownloadProof = onCleanup(@() app.deleteFileIfExists(downloadProof));
            cleanupUrls = onCleanup(@() app.deleteFileIfExists(urlsFile));

            fid = fopen(urlsFile, 'w');
            if fid == -1
                uialert(app.UIFigure, 'Cannot create the temporary URL list.', 'Download Failed');
                updateOutput(app, 'Update download failed: cannot create the temporary URL list.');
                return
            end

            for k = 1:size(selectedRows, 1)
                sceneName = selectedRows{k, 2};
                sceneUrl = app.findSENUpdateSceneUrl(sceneName);
                if isempty(sceneUrl)
                    fclose(fid);
                    uialert(app.UIFigure, ['Missing download URL for scene: ' sceneName], 'Download Failed');
                    updateOutput(app, ['Update download failed: missing URL for scene ' sceneName]);
                    return
                end
                fprintf(fid, '%s\n', sceneUrl);
            end
            fclose(fid);

            app.UpdateSearchStatusLabel.Text = 'Downloading selected images...';
            updateOutput(app, sprintf('Update download started for %d selected Sentinel-1 image(s).', size(selectedRows, 1)));
            drawnow;

            command = app.buildPythonCommand(scriptPath, {'--files', urlsFile, ...
                '--destination', destinationFolder, '--output', outputFile});

            usesDownloadProof = false;
            if ispc
                batchPath = fullfile(projectFolder, 'download_update_images.bat');
                fid = fopen(batchPath, 'w');
                if fid == -1
                    uialert(app.UIFigure, 'Cannot create download_update_images.bat.', 'Download Failed');
                    updateOutput(app, 'Update download failed: cannot create download_update_images.bat.');
                    return
                end
                fprintf(fid, '@echo off \r\n');
                fprintf(fid, 'cd /d "%s"\r\n', currentFolder);
                fprintf(fid, '%s\r\n', command);
                fprintf(fid, 'set DL_STATUS=%%ERRORLEVEL%%\r\n');
                fprintf(fid, 'type nul > "%s"\r\n', downloadProof);
                fprintf(fid, 'exit /b %%DL_STATUS%%\r\n');
                fclose(fid);
                runCommand = sprintf('start "PHASE Update Download" /wait cmd /c ""%s""', batchPath);
                usesDownloadProof = true;
                [status, commandOutput] = system(runCommand);
            elseif isunix && ~ismac
                batchPath = fullfile(projectFolder, 'download_update_images.sh');
                fid = fopen(batchPath, 'w');
                if fid == -1
                    uialert(app.UIFigure, 'Cannot create download_update_images.sh.', 'Download Failed');
                    updateOutput(app, 'Update download failed: cannot create download_update_images.sh.');
                    return
                end
                fprintf(fid, '#!/bin/bash \n');
                fprintf(fid, 'cd "%s"\n', currentFolder);
                fprintf(fid, '%s\n', command);
                fprintf(fid, 'DL_STATUS=$?\n');
                fprintf(fid, 'touch "%s"\n', downloadProof);
                fprintf(fid, 'exit $DL_STATUS\n');
                fclose(fid);
                system(['chmod +x ' app.quoteCommandArgument(batchPath)]);
                usesDownloadProof = true;
                [status, commandOutput] = system(['xterm -e ' app.quoteCommandArgument(batchPath)], '-echo');
            elseif ismac
                batchPath = fullfile(projectFolder, 'download_update_images.sh');
                fid = fopen(batchPath, 'w');
                if fid == -1
                    uialert(app.UIFigure, 'Cannot create download_update_images.sh.', 'Download Failed');
                    updateOutput(app, 'Update download failed: cannot create download_update_images.sh.');
                    return
                end
                fprintf(fid, '#!/bin/bash \n');
                fprintf(fid, 'cd "%s"\n', currentFolder);
                fprintf(fid, '%s\n', command);
                fprintf(fid, 'DL_STATUS=$?\n');
                fprintf(fid, 'touch "%s"\n', downloadProof);
                fprintf(fid, 'exit $DL_STATUS\n');
                fclose(fid);
                system(['chmod +x ' app.quoteCommandArgument(batchPath)]);
                usesDownloadProof = true;
                [status, commandOutput] = system(['open -a Terminal ' app.quoteCommandArgument(batchPath)], '-echo');
            else
                [status, commandOutput] = system(command);
            end

            if usesDownloadProof
                while exist(downloadProof, 'file') == 0
                    pause(1);
                end
                app.deleteFileIfExists(downloadProof);
            end

            if status ~= 0
                app.UpdateSearchStatusLabel.Text = 'Download failed.';
                updateOutput(app, 'Update download failed.');
                uialert(app.UIFigure, commandOutput, 'ASF Download Failed');
                return
            end

            summarySuccessCount = NaN;
            summaryFailedCount = 0;
            summarySkippedCount = 0;
            failedItems = {};

            if exist(outputFile, 'file') == 2
                try
                    downloadSummary = jsondecode(fileread(outputFile));

                    if isfield(downloadSummary, 'successCount')
                        summarySuccessCount = double(downloadSummary.successCount);
                    end

                    if isfield(downloadSummary, 'failedCount')
                        summaryFailedCount = double(downloadSummary.failedCount);
                    end

                    if isfield(downloadSummary, 'skippedCount')
                        summarySkippedCount = double(downloadSummary.skippedCount);
                    end

                    if isfield(downloadSummary, 'failed') && ~isempty(downloadSummary.failed)
                        failedItems = cellstr(string(downloadSummary.failed));
                    end
                catch ME
                    updateOutput(app, ['Could not parse update download summary: ' ME.message]);
                end
            end

            if summaryFailedCount > 0
                app.UpdateSearchStatusLabel.Text = 'Download failed.';
                updateOutput(app, sprintf('Update download failed: %d succeeded, %d skipped, %d failed.', ...
                    max(summarySuccessCount, 0), summarySkippedCount, summaryFailedCount));

                if isempty(failedItems)
                    failureMessage = 'ASF download failed. No files were downloaded.';
                else
                    failureMessage = ['ASF download failed for:' newline strjoin(failedItems, newline)];
                end

                uialert(app.UIFigure, failureMessage, 'ASF Download Failed');
                return
            end

            updateDates = app.getSENRootZipDates(projectFolder);
            if isempty(updateDates)
                app.UpdateSearchStatusLabel.Text = 'Download failed: no new ZIP files found in slaves.';
                updateOutput(app, 'Update download failed: no new Sentinel-1 ZIP files were found in PHASE_Preprocessing/slaves after the download command.');
                uialert(app.UIFigure, ...
                    'Download finished, but no new Sentinel-1 ZIP files were found in PHASE_Preprocessing/slaves. Check the ASF download summary and credentials/network access.', ...
                    'No Downloaded Images');
                return
            end

            if isnan(summarySuccessCount)
                downloadedCount = size(selectedRows, 1);
            else
                downloadedCount = summarySuccessCount;
            end

            app.UpdateSearchStatusLabel.Text = sprintf('Downloaded %d image(s) into slaves.', downloadedCount);
            updateOutput(app, sprintf('Update download completed: %d Sentinel-1 image(s) downloaded into slaves.', downloadedCount));
            success = true;

            if autoRun
                app.runSENUpdateAfterDownload();
            elseif showCompletionDialog
                choice = uiconfirm(app.UIFigure, ...
                    sprintf('Download completed. %d Sentinel-1 image(s) downloaded into slaves.', downloadedCount), ...
                    'Download Completed', ...
                    'Options', {'Close', 'Close and run'}, ...
                    'DefaultOption', 'Close', ...
                    'CancelOption', 'Close', ...
                    'Icon', 'info');
                if strcmp(choice, 'Close and run')
                    app.runSENUpdateAfterDownload();
                end
            end
        end

        function runSENUpdateAfterDownload(app)
            app.TabGroup.SelectedTab = app.RunTab;
            updateOutput(app, 'Starting update preprocessing after download.');
            drawnow;

            try
                app.StartButtonPushed([]);
            catch ME
                updateOutput(app, ['Update preprocessing failed to start: ' ME.message]);
                drawnow;
                uialert(app.UIFigure, ME.message, 'Update Preprocessing Failed');
            end
        end

        function populateSENUpdateImagesTable(app, results)
            if isempty(results)
                app.UpdateImagesTable.Data = {};
                app.UpdateSearchStatusLabel.Text = 'No new matching images found.';
                app.DownloadUpdateImagesButton.Enable = 'off';
                app.DownloadRunUpdateImagesButton.Enable = 'off';
                return
            end

            rows = cell(numel(results), 5);
            for k = 1:numel(results)
                rows{k, 1} = true;
                rows{k, 2} = app.fieldToChar(results(k), 'sceneName');
                rows{k, 3} = app.fieldToChar(results(k), 'startTime');
                rows{k, 4} = app.fieldToChar(results(k), 'platform');
                rows{k, 5} = app.fieldToChar(results(k), 'polarization');
            end

            app.UpdateImagesTable.Data = rows;
            app.UpdateSearchStatusLabel.Text = sprintf('Found %d new matching image(s).', size(rows, 1));
            app.DownloadUpdateImagesButton.Enable = 'on';
            app.DownloadRunUpdateImagesButton.Enable = 'on';
        end

        function sceneUrl = findSENUpdateSceneUrl(app, sceneName)
            sceneUrl = '';
            results = app.update_search_results_SEN;
            if isempty(results)
                return
            end

            for k = 1:numel(results)
                currentScene = app.fieldToChar(results(k), 'sceneName');
                if strcmp(currentScene, sceneName)
                    sceneUrl = app.fieldToChar(results(k), 'url');
                    if strcmp(sceneUrl, '-')
                        sceneUrl = '';
                    end
                    return
                end
            end
        end

        function setAllSENUpdateImagesSelected(app, selected)
            tableData = app.UpdateImagesTable.Data;
            if isempty(tableData)
                return
            end

            for k = 1:size(tableData, 1)
                tableData{k, 1} = selected;
            end
            app.UpdateImagesTable.Data = tableData;
        end

        function parametersText = formatSENUpdateParameters(app, params)
            parametersText = {'File type: L1 Single Look Complex (SLC)'; ...
                'Beam mode: IW'; ...
                'Polarization: empty (all)'; ...
                ['Direction: ' app.fieldToChar(params, 'flightDirection')]; ...
                'Subtype: All (Sentinel-1A/B/C/D)'; ...
                ['Path Start/End: ' app.fieldToChar(params, 'pathNumber')]; ...
                ['Frame Start/End: ' app.fieldToChar(params, 'frameNumber')]; ...
                ['Reference: ' app.fieldToChar(params, 'referenceScene')]};
        end

        function command = buildPythonCommand(app, scriptPath, arguments)
            parts = [{app.python_SEN}, {scriptPath}, arguments];
            quotedParts = cell(size(parts));
            for k = 1:numel(parts)
                quotedParts{k} = app.quoteCommandArgument(parts{k});
            end
            command = strjoin(quotedParts, ' ');
        end

        function quoted = quoteCommandArgument(app, value)
            value = char(string(value));
            value = strrep(value, '"', '""');
            quoted = ['"' value '"'];
        end

        function textValue = fieldToChar(app, data, fieldName)
            textValue = '-';
            if isstruct(data) && isfield(data, fieldName) && ~isempty(data.(fieldName))
                value = data.(fieldName);
                if isnumeric(value)
                    textValue = num2str(value);
                elseif isstring(value)
                    textValue = char(value);
                elseif ischar(value)
                    textValue = value;
                elseif iscell(value)
                    textValue = strjoin(cellfun(@char, value, 'UniformOutput', false), ', ');
                else
                    textValue = char(string(value));
                end
            end
        end

        function deleteFileIfExists(app, filePath)
            if exist(filePath, 'file') == 2
                delete(filePath);
            end
        end

        % check if logged in
        function loggedIn = isEarthdataLoggedIn(app)

            loggedIn = strcmp(app.SignedInLabel.Visible, "on");

        end

        function restoreLoginState(app)

            appPath = phase_preprocessing_beta.projectRoot();
            resultFile = fullfile(appPath, "downloadasf", "login_result.json");
            requestFile = fullfile(appPath, "downloadasf", "login_request.json");

            if ~exist(resultFile, "file") || ~exist(requestFile, "file")
                return;
            end

            result = jsondecode(fileread(resultFile));

            if ~isfield(result, "status") || ~strcmp(string(result.status), "success")
                return;
            end

            request = jsondecode(fileread(requestFile));

            if ~isfield(request, "username") || strlength(string(request.username)) == 0
                return;
            end

            username = string(request.username);

            app.EarthdataUsernameEditField.Visible = "off";
            app.EarthdataPasswordEditField.Visible = "off";
            app.LoginButton.Visible = "off";

            app.SignedInLabel.Text = "Signed in as " + username;
            app.SignedInLabel.Visible = "on";

            app.SignOutButton.Visible = "on";
            app.LoginFeedbackLabel.Text = "";

        end

                function initializeCSKMap(app)
            % 1. Pulisce completamente il pannello e rigenera gli assi geometrici
            delete(app.Panel.Children);
            gx = geoaxes(app.Panel);
            geobasemap(gx, 'satellite');
            gx.Interactions = [panInteraction; zoomInteraction];

            % 2. Controlla lo stato attuale della cartella slaves per Cosmo-SkyMed
            currentFolder = phase_preprocessing_beta.projectRoot();
            project_path_full = fullfile(currentFolder, 'PHASE_Preprocessing');
            csk_files = dir(fullfile(project_path_full, 'slaves', '**', '*.h5'));

            default_pos = [];

            % 3. Se sono presenti file (anche se scaricati dopo l'avvio), gestisce il primo
            if ~isempty(csk_files)
                curr_file = fullfile(csk_files(1).folder, csk_files(1).name);
                try
                    % Estrazione degli attributi di georeferenziazione dall'HDF5
                    tl = h5readatt(curr_file, '/', 'Estimated Top Left Geodetic Coordinates');
                    tr = h5readatt(curr_file, '/', 'Estimated Top Right Geodetic Coordinates');
                    br = h5readatt(curr_file, '/', 'Estimated Bottom Right Geodetic Coordinates');
                    bl = h5readatt(curr_file, '/', 'Estimated Bottom Left Geodetic Coordinates');

                    lats = [tl(1), tr(1), br(1), bl(1)];
                    lons = [tl(2), tr(2), br(2), bl(2)];

                    % Disegna il footprint blu del primo file trovato
                    hold(gx, 'on');
                    geoplot(gx, [lats lats(1)], [lons lons(1)], 'Color', 'b', 'LineWidth', 1.5);
                    hold(gx, 'off');

                    % Imposta i limiti della mappa attorno al footprint
                    geolimits(gx, [min(lats)-0.05, max(lats)+0.05], [min(lons)-0.05, max(lons)+0.05]);

                    % Calcola il box rosso iniziale centrato all'80% del footprint
                    lon_min_fp = min(lons); lon_max_fp = max(lons);
                    lat_min_fp = min(lats); lat_max_fp = max(lats);

                    w = (lon_max_fp - lon_min_fp) * 0.8;
                    h = (lat_max_fp - lat_min_fp) * 0.8;
                    start_lon = lon_min_fp + (lon_max_fp - lon_min_fp) * 0.1;
                    start_lat = lat_min_fp + (lat_max_fp - lat_min_fp) * 0.1;

                    default_pos = [start_lat, start_lon, h, w];
                catch
                    disp(['Impossibile leggere i metadati HDF5 da: ', csk_files(1).name]);
                end
            end

            % 4. Disegna il rettangolo interattivo dell'AOI (regione di interesse)
            if ~isempty(default_pos)
                app.roi_CSK = drawrectangle(gx, 'Position', default_pos, 'Color', 'r', 'FaceAlpha', 0.2);
                app.updateCSKCoords(app.roi_CSK); % Sincronizza i campi numerici della UI
            else
                % Se non ci sono ancora file, genera un ROI editabile standard (es. sull'Italia)
                app.roi_CSK = drawrectangle(gx, 'Position', [42, 12, 1, 1], 'Color', 'r', 'FaceAlpha', 0.2);
                app.updateCSKCoords(app.roi_CSK);
            end

            % Associa il listener per aggiornare i parametri numerici durante il trascinamento
            addlistener(app.roi_CSK, 'ROIMoved', @(src, event) app.updateCSKCoords(src));
        end


        function updateCSKCoords(app, roi_source)
            if isvalid(roi_source)
                % FIXED: GeoAxes outputs [LatMin, LonMin, LatHeight, LonWidth]
                pos = roi_source.Position;

                latMin = pos(1);
                lonMin = pos(2);
                latMax = pos(1) + pos(3);
                lonMax = pos(2) + pos(4);

                % Update UI Fields immediately
                app.MinlongitudeEditField_2.Value = lonMin;
                app.MinlatitudeEditField_2.Value = latMin;
                app.MaxlongitudeEditField_2.Value = lonMax;
                app.MaxlatitudeEditField_2.Value = latMax;

                % Update internal background variables
                app.lon_min_CSK = lonMin;
                app.lat_min_CSK = latMin;
                app.lon_max_CSK = lonMax;
                app.lat_max_CSK = latMax;
            end
        end


                function initializeSENMap(app)
            % 1. Pulisce completamente il pannello e rigenera gli assi
            delete(app.SENMapPanel.Children);
            gx = geoaxes(app.SENMapPanel);
            geobasemap(gx, 'satellite');
            gx.Interactions = [panInteraction; zoomInteraction];

            % 2. Controlla lo stato attuale della cartella slaves
            currentFolder = phase_preprocessing_beta.projectRoot();
            project_path_full = fullfile(currentFolder, 'PHASE_Preprocessing');
            sen_files = dir(fullfile(project_path_full, 'slaves', '**', '*.zip'));

            default_pos = [];

            % 3. Se sono apparsi file (anche dopo il lancio dell'app), legge il primo
            if ~isempty(sen_files)
                curr_file = fullfile(sen_files(1).folder, sen_files(1).name);
                try
                    % --- FAST JAVA ZIP TRICK ---
                    zipObj = java.util.zip.ZipFile(curr_file);
                    entries = zipObj.entries();
                    manifest_text = '';
                    while entries.hasMoreElements()
                        entry = entries.nextElement();
                        if endsWith(char(entry.getName()), 'manifest.safe')
                            stream = zipObj.getInputStream(entry);
                            scanner = java.util.Scanner(stream).useDelimiter('\A');
                            if scanner.hasNext()
                                manifest_text = char(scanner.next());
                            end
                            scanner.close();
                            break;
                        end
                    end
                    zipObj.close();

                    % PARSE XML COORDINATES
                    coords_tokens = regexp(manifest_text, '<gml:coordinates>([^<]+)</gml:coordinates>', 'tokens', 'once');
                    if ~isempty(coords_tokens)
                        coords_str = strtrim(coords_tokens{1});
                        pairs = split(coords_str, ' ');
                        lats = zeros(1, length(pairs));
                        lons = zeros(1, length(pairs));
                        for k = 1:length(pairs)
                            latlon = split(pairs{k}, ',');
                            lats(k) = str2double(latlon{1});
                            lons(k) = str2double(latlon{2});
                        end

                        % Disegna il footprint blu del primo file trovato
                        hold(gx, 'on');
                        geoplot(gx, [lats lats(1)], [lons lons(1)], 'Color', 'b', 'LineWidth', 1.5);
                        hold(gx, 'off');

                        geolimits(gx, [min(lats)-0.05, max(lats)+0.05], [min(lons)-0.05, max(lons)+0.05]);

                        % Calcola il box rosso centrato all'80%
                        lon_min_fp = min(lons); lon_max_fp = max(lons);
                        lat_min_fp = min(lats); lat_max_fp = max(lats);
                        w = (lon_max_fp - lon_min_fp) * 0.8;
                        h = (lat_max_fp - lat_min_fp) * 0.8;
                        start_lon = lon_min_fp + (lon_max_fp - lon_min_fp) * 0.1;
                        start_lat = lat_min_fp + (lat_max_fp - lat_min_fp) * 0.1;
                        default_pos = [start_lat, start_lon, h, w];
                    end
                catch
                    disp(['Impossibile leggere i metadati da: ', sen_files(1).name]);
                end
            end

            % 4. Disegna il rettangolo interattivo dell'AOI
            if ~isempty(default_pos)
                app.roi_SEN = drawrectangle(gx, 'Position', default_pos, 'Color', 'r', 'FaceAlpha', 0.2);
                app.updateSENCoords(app.roi_SEN);
            else
                % Se la cartella è ancora vuota, crea un rettangolo di default (es. sull'Italia)
                % per non lasciare l'interfaccia senza un ROI attivo modificabile
                app.roi_SEN = drawrectangle(gx, 'Position', [42, 12, 1, 1], 'Color', 'r', 'FaceAlpha', 0.2);
                app.updateSENCoords(app.roi_SEN);
            end

            % Ricollega il listener per aggiornare le coordinate numeriche nella UI
            addlistener(app.roi_SEN, 'ROIMoved', @(src, event) app.updateSENCoords(src));
        end


        function updateSENCoords(app, roi_source)
            if isvalid(roi_source)
                % GeoAxes outputs [LatMin, LonMin, LatHeight, LonWidth]
                pos = roi_source.Position;

                latMin = pos(1);
                lonMin = pos(2);
                latMax = pos(1) + pos(3);
                lonMax = pos(2) + pos(4);

                % Update Sentinel-1 UI Fields immediately
                app.MinlongitudeEditField.Value = lonMin;
                app.MinlatitudeEditField.Value = latMin;
                app.MaxlongitudeEditField.Value = lonMax;
                app.MaxlatitudeEditField.Value = latMax;

                % Update Sentinel-1 internal variables
                app.lon_min_SEN = lonMin;
                app.lat_min_SEN = latMin;
                app.lon_max_SEN = lonMax;
                app.lat_max_SEN = latMax;
            end
        end


         function refreshImportedImagesTable(app)
              % Refresh the table with the .h5 files currently in slaves.
              currentFolder = phase_preprocessing_beta.projectRoot();
              slavesFolder = fullfile(currentFolder, 'PHASE_Preprocessing', 'slaves');
              h5list = dir(fullfile(slavesFolder, '**', '*.h5'));
              if isempty(h5list)
                  app.ImportedImagesTable.Data = {};
                  return
              end
              rows = cell(numel(h5list), 4);
              for k = 1:numel(h5list)
                  [dstr, tipo] = parseCSKInfo(app, h5list(k).name);
                  rows{k,1} = h5list(k).name;
                  rows{k,2} = dstr;
                  rows{k,3} = tipo;
                  rows{k,4} = 'OK';
              end
              app.ImportedImagesTable.Data = rows;
              app.ImportedImagesTable.ColumnName = {'File','Date','Type','Status'};
          end

              function [dstr, tipo] = parseCSKInfo(app, fname)
              % Detect generation (CSK/CSG) and acquisition date from the file name.
              if startsWith(upper(fname), 'CSG')
                  tipo = 'CSG';
              else
                  tipo = 'CSK';
              end
              tok = regexp(fname, '(\d{8})', 'tokens', 'once');
              if isempty(tok)
                  dstr = '-';
              else
                  d = tok{1};
                  dstr = sprintf('%s-%s-%s', d(1:4), d(5:6), d(7:8));
              end
          end

        function refreshImportedSENImagesTable(app)
            if isempty(app.ImportedSENImagesTable) || ~isvalid(app.ImportedSENImagesTable)
                return;
            end

            currentFolder = phase_preprocessing_beta.projectRoot();
            projectFolder = fullfile(currentFolder, 'PHASE_Preprocessing');
            slavesFolder = fullfile(projectFolder, 'slaves');

            if exist(slavesFolder, 'dir') ~= 7
                app.ImportedSENImagesTable.Data = {};
                return;
            end

            files = dir(fullfile(slavesFolder, '*.zip'));

            if isempty(files)
                app.ImportedSENImagesTable.Data = {};
                return;
            end

            data = cell(numel(files), 3);

            for k = 1:numel(files)
                fileName = files(k).name;
                dateText = app.parseSENDateFromFilename(fileName);

                if isempty(dateText)
                    statusText = 'Ready';
                elseif app.hasSENProcessedSlaveDate(projectFolder, dateText)
                    statusText = 'Already processed';
                else
                    statusText = 'Ready';
                end

                data{k, 1} = fileName;
                data{k, 2} = dateText;
                data{k, 3} = statusText;
            end

            app.ImportedSENImagesTable.Data = data;
        end
    end


    % Callbacks that handle component events
    methods (Access = public)

        % Code that executes after component creation
        function startupFcn(app)

            % Set the working directory to the PHASE project root
            currentFolder = phase_preprocessing_beta.projectRoot();
            cd(currentFolder);

            % Ensure Mapping Toolbox and Image Processing Toolbox are installed
            installedToolboxes = {ver().Name};
            if ~any(strcmp(installedToolboxes, 'Mapping Toolbox')) || ~any(strcmp(installedToolboxes, 'Image Processing Toolbox'))
                msg = 'Missing required toolboxes! Please install Mapping Toolbox and Image Processing Toolbox before running this app.';
                uialert(app.UIFigure, msg, 'Missing Dependencies');

                % Display message in the log windows just in case
                app.MessagesTextArea.Value = {'ERROR: Missing Mapping Toolbox or Image Processing Toolbox.'};
                app.MessagesTextArea_2.Value = {'ERROR: Missing Mapping Toolbox or Image Processing Toolbox.'};
                return; % Exit the startup function early
            end

            % Initially hide both the label and the edit field
            app.PythonEnvironmentLabel.Visible = 'off';
            app.CustomPythonEnvironmentEditField.Visible = 'off';
            app.CustomPythonEnvironmentLabel.Visible = 'off';
            app.CustomPythonEnvironmentEditField_2.Visible = 'off';

            % AUTO-CONFIG (PHASE installer): default panel SEN + path Python/GPT
            try
                % 1. Default Sentinel1 panel visible
                if isprop(app, 'ConstellationSwitch') && isvalid(app.ConstellationSwitch)
                    app.ConstellationSwitch.Value = 'Sentinel1';
                end
                if isprop(app, 'Sentinel1Panel') && isvalid(app.Sentinel1Panel)
                    app.Sentinel1Panel.Visible = 'on';
                end
                if isprop(app, 'CosmoSkyMedPanel') && isvalid(app.CosmoSkyMedPanel)
                    app.CosmoSkyMedPanel.Visible = 'off';
                end
                app.constellation = 'SEN';
                drawnow;

                % 2. Set diretto dei 2 path da input_preprocessing.mat (Python+GPT)
                if exist('./PHASE_Preprocessing/input_preprocessing.mat', 'file') == 2
                    cfg = load('./PHASE_Preprocessing/input_preprocessing.mat');
                    if isfield(cfg, 'python')
                        % SEN side
                        if isprop(app, 'CustomPythonEnvironmentEditField') && isvalid(app.CustomPythonEnvironmentEditField)
                            app.CustomPythonEnvironmentEditField.Value = cfg.python;
                            app.CustomPythonEnvironmentEditField.Visible = 'on';
                        end
                        if isprop(app, 'PythonEnvironmentDropDown') && isvalid(app.PythonEnvironmentDropDown)
                            app.PythonEnvironmentDropDown.Value = 'Other';
                        end
                        if isprop(app, 'PythonEnvironmentLabel') && isvalid(app.PythonEnvironmentLabel)
                            app.PythonEnvironmentLabel.Visible = 'on';
                        end
                        app.python_SEN = cfg.python;
                        % CSK side (stesso python_exe)
                        if isprop(app, 'CustomPythonEnvironmentEditField_2') && isvalid(app.CustomPythonEnvironmentEditField_2)
                            app.CustomPythonEnvironmentEditField_2.Value = cfg.python;
                            app.CustomPythonEnvironmentEditField_2.Visible = 'on';
                        end
                        if isprop(app, 'PythonEnvironmentDropDown_2') && isvalid(app.PythonEnvironmentDropDown_2)
                            app.PythonEnvironmentDropDown_2.Value = 'Other';
                        end
                        app.python_CSK = cfg.python;
                    end
                    if isfield(cfg, 'gptbin_path')
                        if isprop(app, 'PathEditField') && isvalid(app.PathEditField)
                            app.PathEditField.Value = cfg.gptbin_path;
                        end
                        app.gptbin_path_SEN = cfg.gptbin_path;
                        if isprop(app, 'PathEditField_2') && isvalid(app.PathEditField_2)
                            app.PathEditField_2.Value = cfg.gptbin_path;
                        end
                        app.gptbin_path_CSK = cfg.gptbin_path;
                    end
                end
            catch
                % non blocca lo startup se qualcosa fallisce
            end

            % Render the CSK map
            drawnow;
            app.initializeCSKMap();
            app.initializeSENMap();



            % Initialize Downloader tab
            initializeDownloaderMap(app);
            % Default values for params
            app.FilterPanel.Visible = "off";
            app.SignOutButton.Visible = "off";
            app.SignedInLabel.Visible = "off";
            app.initializeParams();
            app.PreviewPanel.Visible = "off";
            app.DownloadSelectedButton.Enable = "off";
            app.PreviewSelectedLabel.Text = "Selected: 0";
            app.DownloadProgressGauge.Visible = "off";
            app.DownloadProgressLabel.Visible = "off";
            app.DownloadProgressGauge.Value = 0;
            app.DownloadProgressLabel.Text = "";
            app.restoreLoginState();
            refreshImportedImagesTable(app);
            app.refreshImportedSENImagesTable();
            app.UpdateSearchEndDateDatePicker.Value = datetime('today');
            app.setSENUpdateControlsVisible(false);

        end

        % Value changed function: ConstellationSwitch
        function ConstellationSwitchValueChanged(app, event)
           selectedOption = app.ConstellationSwitch.Value;

            if selectedOption == "Sentinel1"
                app.Sentinel1Panel.Visible = 'on';
                app.CosmoSkyMedPanel.Visible = 'off';
                app.constellation = 'SEN';
            elseif selectedOption == "CosmoSkyMed"
                app.Sentinel1Panel.Visible = 'off';
                app.CosmoSkyMedPanel.Visible = 'on';
                app.constellation = 'CSK';
            end

        end

        % Callback function
        function RunningOperativeSystemButtonGroupSelectionChanged(app, event)
            selectedButtonText = app.RunningOperativeSystemButtonGroup.SelectedObject.Text;

            if strcmp(selectedButtonText, "Windows")
                app.os_SEN = 0;
            elseif strcmp(selectedButtonText, "Linux")
                app.os_SEN = 1;
            end

        end

        % Value changed function: CustomPythonEnvironmentEditField
        function CustomPythonEnvironmentEditFieldValueChanged(app, event)
            value = app.CustomPythonEnvironmentEditField.Value;

            app.python_SEN = value;

            selectedValue = app.PythonEnvironmentDropDown.Value;

            if strcmp(selectedValue, 'Other')
                % If "other" is selected, use the value from the EditField
                app.python_SEN = app.CustomPythonEnvironmentEditField.Value;
            else
                % Use the selected value from the dropdown menu
                app.python_SEN = selectedValue;
            end

        end

        % Button pushed function: SaveButton, SaveButton_2
        function SaveButtonPushed(app, event)

            % Specify the filename for saving
            constellation = app.constellation;
            filename = './PHASE_Preprocessing/input_preprocessing.mat';  % specify the filename for saving

            if app.constellation == "SEN"

                python = app.python_SEN;
                master_date = app.master_date_SEN;
                auto_master = app.auto_master_SEN;
                master_processing = app.master_processing_SEN;
                polarisation = app.polarisation_SEN;
                lon_min = app.lon_min_SEN;
                lat_min = app.lat_min_SEN;
                lon_max = app.lon_max_SEN;
                lat_max = app.lat_max_SEN;
                slaves_removal = app.slaves_removal_SEN;
                dem_name = app.dem_name_SEN;
                dem_file = app.dem_file_SEN;
                dem_name_coreg = app.dem_name_coreg_SEN;
                dem_file_coreg = app.dem_file_coreg_SEN;
                dem_resampling = app.dem_resampling_SEN;
                first_step = app.first_step_SEN;
                coherence_tc = app.coherence_tc_SEN;
                epsg_code = app.epsg_code_SEN;
                gptbin_path = app.gptbin_path_SEN;
                cpu = app.cpu_SEN;
                cache = app.cache_SEN;

                try

                    % Delete the existing .mat file, if it exists
                    if exist(filename, 'file') == 2
                        delete(filename);
                    end

                    save(filename, 'constellation', 'python', ...
                        'master_date', 'auto_master', 'master_processing', ...
                        'polarisation', 'lon_min', ...
                        'lat_min', 'lon_max', 'lat_max', 'slaves_removal', ...
                        'dem_name', 'dem_file', ...
                        'dem_name_coreg', 'dem_file_coreg', 'dem_resampling', 'first_step', 'coherence_tc', ...
                        'epsg_code', 'gptbin_path', 'cpu', 'cache', '-mat');

                    % If saving is successful, set the lamp to green
                    app.StatusLamp.Color = [0, 1, 0]; % Green color

                catch exception

                    % If saving fails, set the lamp to red
                    app.StatusLamp.Color = [1, 0, 0]; % Red color
                    % Display the error message
                    updateOutput(app, ['Error saving file: ' exception.message]);

                end

            elseif app.constellation == "CSK"

                python = app.python_CSK;
                master_date = app.master_date_CSK;
                auto_master = app.auto_master_CSK;
                master_processing = app.master_processing_CSK;
                lon_min = app.lon_min_CSK;
                lat_min = app.lat_min_CSK;
                lon_max = app.lon_max_CSK;
                lat_max = app.lat_max_CSK;
                slaves_removal = app.slaves_removal_CSK;
                dem_name = app.dem_name_CSK;
                dem_file = app.dem_file_CSK;
                first_step = app.first_step_CSK;
                num_gcp = app.num_gcp_CSK;
                coherence_tc = app.coherence_tc_CSK;
                epsg_code = app.epsg_code_CSK;
                gptbin_path = app.gptbin_path_CSK;
                cpu = app.cpu_CSK;
                cache = app.cache_CSK;

                try

                    % Delete the existing .mat file, if it exists
                    if exist(filename, 'file') == 2
                        delete(filename);
                    end

                    save(filename, 'constellation', 'python', ...
                        'master_date', 'auto_master', 'master_processing', 'lon_min', ...
                        'lat_min', 'lon_max', 'lat_max', 'slaves_removal', 'dem_name', 'dem_file', ...
                        'first_step', 'num_gcp', 'coherence_tc', 'epsg_code', 'gptbin_path', 'cpu', 'cache', '-mat');

                    % If saving is successful, set the lamp to green
                    app.StatusLamp_2.Color = [0, 1, 0]; % Green color

                catch exception

                    % If saving fails, set the lamp to red
                    app.StatusLamp_2.Color = [1, 0, 0]; % Red color
                    % Display the error message
                    updateOutput(app, ['Error saving file: ' exception.message]);

                end

            end


        end

        % Value changed function: MasterprocessingCheckBox
        function MasterprocessingCheckBoxValueChanged(app, event)
            value = app.MasterprocessingCheckBox.Value;

            if value
                app.master_processing_SEN = 0;
            else
                app.master_processing_SEN = 1;
            end

        end

        % Value changed function: PolarisationDropDown
        function PolarisationDropDownValueChanged(app, event)
            value = app.PolarisationDropDown.Value;

            app.polarisation_SEN = value;

        end

        % Value changed function: MinlongitudeEditField
        function MinlongitudeEditFieldValueChanged(app, event)
            value = app.MinlongitudeEditField.Value;

            app.lon_min_SEN = value;

        end

        % Value changed function: MaxlongitudeEditField
        function MaxlongitudeEditFieldValueChanged(app, event)
            value = app.MaxlongitudeEditField.Value;

            app.lon_max_SEN = value;

        end

        % Value changed function: MaxlatitudeEditField
        function MaxlatitudeEditFieldValueChanged(app, event)
            value = app.MaxlatitudeEditField.Value;

            app.lat_max_SEN = value;

        end

        % Value changed function: MinlatitudeEditField
        function MinlatitudeEditFieldValueChanged(app, event)
            value = app.MinlatitudeEditField.Value;

            app.lat_min_SEN = value;

        end

        % Value changed function: MasterdateDatePicker
        function MasterdateDatePickerValueChanged(app, event)
            value = app.MasterdateDatePicker.Value;

            date_txt = datetime(value, 'Format', 'ddMMyyyy');
            app.master_date_SEN = date_txt;

        end

        % Value changed function: FirststepDropDown
        function FirststepDropDownValueChanged(app, event)
            value = app.FirststepDropDown.Value;

            app.first_step_SEN = value;

        end

        % Value changed function: DEMresamplingmethodDropDown
        function DEMresamplingmethodDropDownValueChanged(app, event)
            value = app.DEMresamplingmethodDropDown.Value;

            app.dem_resampling_SEN = value;

        end

        % Value changed function: DEMcoregpathEditField
        function DEMcoregpathEditFieldValueChanged(app, event)
            value = app.DEMcoregpathEditField.Value;

            app.dem_file_coreg_SEN = value;

        end

        % Value changed function: DEMifgpathEditField
        function DEMifgpathEditFieldValueChanged(app, event)
            value = app.DEMifgpathEditField.Value;

            app.dem_file_SEN = value;

        end

        % Value changed function: DEMcoregistrationDropDown
        function DEMcoregistrationDropDownValueChanged(app, event)
            value = app.DEMcoregistrationDropDown.Value;

            app.dem_name_coreg_SEN = value;

        end

        % Value changed function: DEMinterferogramDropDown
        function DEMinterferogramDropDownValueChanged(app, event)
            value = app.DEMinterferogramDropDown.Value;

            app.dem_name_SEN = value;

        end

        % Value changed function: SlavesremovalafterprocessingCheckBox
        function SlavesremovalafterprocessingCheckBoxValueChanged(app, event)
            value = app.SlavesremovalafterprocessingCheckBox.Value;

            if value
                app.slaves_removal_SEN = 0;
            else
                app.slaves_removal_SEN = 1;
            end

        end

        % Value changed function: EPSGcodeEditField
        function EPSGcodeEditFieldValueChanged(app, event)
            value = app.EPSGcodeEditField.Value;

            app.epsg_code_SEN = value;

        end

        % Value changed function: TerraincorrectedCoherenceandLIACheckBox
        function TerraincorrectedCoherenceandLIACheckBoxValueChanged(app, event)
            value = app.TerraincorrectedCoherenceandLIACheckBox.Value;

            if value
                app.coherence_tc_SEN = 0;
            else
                app.coherence_tc_SEN = 1;
            end

        end

        % Value changed function: CacheEditField
        function CacheEditFieldValueChanged(app, event)
            value = app.CacheEditField.Value;

            app.cache_SEN = value;

        end

        % Value changed function: CPUEditField
        function CPUEditFieldValueChanged(app, event)
            value = app.CPUEditField.Value;

            app.cpu_SEN = value;

        end

        % Value changed function: PathEditField
        function PathEditFieldValueChanged(app, event)
            value = app.PathEditField.Value;

            app.gptbin_path_SEN = value;

        end

        % Callback function
        function RunningOperativeSystemButtonGroup_2SelectionChanged(app, event)
            selectedButtonText = app.RunningOperativeSystemButtonGroup_2.SelectedObject.Text;

            if strcmp(selectedButtonText, "Windows")
                app.os_CSK = 0;
            elseif strcmp(selectedButtonText, "Linux")
                app.os_CSK = 1;
            end

        end

        % Value changed function: CustomPythonEnvironmentEditField_2
        function CustomPythonEnvironmentEditField_2ValueChanged(app, event)
            value = app.CustomPythonEnvironmentEditField_2.Value;

            app.python_CSK = value;

            selectedValue = app.PythonEnvironmentDropDown_2.Value;

            if strcmp(selectedValue, 'Other')
                % If "other" is selected, use the value from the EditField
                app.python_CSK = app.CustomPythonEnvironmentEditField_2.Value;
            else
                % Use the selected value from the dropdown menu
                app.python_CSK = selectedValue;
            end

        end

        % Value changed function: MinlatitudeEditField_2
        function MinlatitudeEditField_2ValueChanged(app, event)
            value = app.MinlatitudeEditField_2.Value;

            app.lat_min_CSK = value;

        end

        % Value changed function: MaxlatitudeEditField_2
        function MaxlatitudeEditField_2ValueChanged(app, event)
            value = app.MaxlatitudeEditField_2.Value;

            app.lat_max_CSK = value;

        end

        % Value changed function: MaxlongitudeEditField_2
        function MaxlongitudeEditField_2ValueChanged(app, event)
            value = app.MaxlongitudeEditField_2.Value;

            app.lon_max_CSK = value;

        end

        % Value changed function: MinlongitudeEditField_2
        function MinlongitudeEditField_2ValueChanged(app, event)
            value = app.MinlongitudeEditField_2.Value;

            app.lon_min_CSK = value;

        end

        % Value changed function: MasterprocessingCheckBox_2
        function MasterprocessingCheckBox_2ValueChanged(app, event)
            value = app.MasterprocessingCheckBox_2.Value;

            if value
                app.master_processing_CSK = 0;
            else
                app.master_processing_CSK = 1;
            end

        end

        % Value changed function: MasterdateDatePicker_2
        function MasterdateDatePicker_2ValueChanged(app, event)
            value = app.MasterdateDatePicker_2.Value;

            date_txt = datetime(value, 'Format', 'ddMMyyyy');
            app.master_date_CSK = date_txt;

        end

        % Value changed function: CoregistrationGCPsnumberEditField
        function CoregistrationGCPsnumberEditFieldValueChanged(app, event)
            value = app.CoregistrationGCPsnumberEditField.Value;

            app.num_gcp_CSK = value;

        end

        % Value changed function: FirststepDropDown_2
        function FirststepDropDown_2ValueChanged(app, event)
            value = app.FirststepDropDown_2.Value;

            app.first_step_CSK = value;

        end

        % Value changed function: DEMifgpathEditField_2
        function DEMifgpathEditField_2ValueChanged(app, event)
            value = app.DEMifgpathEditField_2.Value;

            app.dem_file_CSK = value;

        end

        % Value changed function: DEMinterferogramDropDown_2
        function DEMinterferogramDropDown_2ValueChanged(app, event)
            value = app.DEMinterferogramDropDown_2.Value;

            app.dem_name_CSK = value;

        end

        % Value changed function: SlavesremovalafterprocessingCheckBox_2
        function SlavesremovalafterprocessingCheckBox_2ValueChanged(app, event)
            value = app.SlavesremovalafterprocessingCheckBox_2.Value;

            if value
                app.slaves_removal_CSK = 0;
            else
                app.slaves_removal_CSK = 1;
            end

        end

        % Value changed function: TerraincorrectedCoherenceandLIACheckBox_2
        function TerraincorrectedCoherenceandLIACheckBox_2ValueChanged(app, event)
            value = app.TerraincorrectedCoherenceandLIACheckBox_2.Value;

            if value
                app.coherence_tc_CSK = 0;
            else
                app.coherence_tc_CSK = 1;
            end

        end

        % Value changed function: EPSGcodeEditField_2
        function EPSGcodeEditField_2ValueChanged(app, event)
            value = app.EPSGcodeEditField_2.Value;

            app.epsg_code_CSK = value;

        end

        % Value changed function: CacheEditField_2
        function CacheEditField_2ValueChanged(app, event)
            value = app.CacheEditField_2.Value;

            app.cache_CSK = value;

        end

        % Value changed function: CPUEditField_2
        function CPUEditField_2ValueChanged(app, event)
            value = app.CPUEditField_2.Value;

            app.cpu_CSK = value;

        end

        % Value changed function: PathEditField_2
        function PathEditField_2ValueChanged(app, event)
            value = app.PathEditField_2.Value;

            app.gptbin_path_CSK = value;

        end

        % Button pushed function: StartButton, StartButton_2
        function StartButtonPushed(app, event)

            app.MessagesTextArea.Value = '';
            app.MessagesTextArea_2.Value = '';
            SaveButtonPushed(app,event);

            if app.constellation == "SEN"

                condition = 0;

                while condition == 0

                    % Display a message in the UIAxes
                    updateOutput(app, 'Button pressed! Running the script...');

                    % Change lamp color to yellow (running)
                    app.PreprocessingstatusLamp.Color = [1, 1, 0]; % Yellow

                    try

                        %% SENTINEL-1 PREPROCESSING CODE ----------------------

                        pwd;
                        prep_folder = pwd;
                        addpath(genpath(prep_folder));

                        %% ------------------ IMPORT OF THE REQUIRED VARIABLES --------------------


                        % DEFINE THE TYPE OF SLASH
                        par = filesep;

                        % READ OF THE INPUT VARIABLES FROM .MAT FILE
                        load(strcat('.', par, 'PHASE_Preprocessing', par, 'input_preprocessing.mat'), 'python', ...
                                                'master_date', 'master_processing', 'auto_master', ...
                                                'polarisation', 'lon_min', 'lat_min', 'lon_max', 'lat_max', 'slaves_removal', 'dem_name', ...
                                                'dem_file', 'dem_name_coreg', 'dem_file_coreg', 'dem_resampling', 'first_step', 'coherence_tc', ...
                                                'epsg_code', 'gptbin_path', 'cpu', 'cache');

                        % PROJECT FOLDER
                        project_path_full = strcat(prep_folder, par, 'PHASE_Preprocessing');

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
                            graphsfolder = strcat('GRAPHSFOLDER=', project_path_full, par, 'snap2stamps', par, 'graphs');

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
                            step_selector = ('SEN_master_selector.py project_master.conf');
                            step_selector_cmd = [python space step_selector];
                            step1_master = ('SEN_splitting_master.py project_master.conf');
                            step_master_2 = [python space step1_master];

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
                                fprintf(f_snap2stamps_master,'#!/bin/bash \n');
                                fprintf(f_snap2stamps_master,'cd "%s"\n',path_1_master);
                                fprintf(f_snap2stamps_master,'%s\n',step_selector_cmd);
                                fprintf(f_snap2stamps_master,'%s\n',step_master_2);
                                fclose(f_snap2stamps_master);

                                path_2_master = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'));
                                chmod = ['chmod +x' space path_2_master];
                                system(chmod);
                                phase_preprocessing_beta.runCommandLive(app, path_2_master, ...
                                    'Master selection and preparation', 3, 18, 1, 1);

                            elseif ismac
                                f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                j = 0;
                                while f_snap2stamps_master == -1 && j < 20
                                    f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                    j = j+1;
                                end
                                fprintf(f_snap2stamps_master,'#!/bin/bash \n');
                                fprintf(f_snap2stamps_master,'cd "%s"\n',path_1_master);
                                fprintf(f_snap2stamps_master,'%s\n',step_selector_cmd);
                                fprintf(f_snap2stamps_master,'%s\n',step_master_2);
                                fclose(f_snap2stamps_master);

                                path_2_master = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'));
                                chmod = ['chmod +x' space path_2_master];
                                system(chmod);
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
                        graphsfolder = strcat('GRAPHSFOLDER=', project_path_full, par, 'snap2stamps', par, 'graphs'); % path of the graphs folder of snap2stamps
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
                        dp = ('::');
                        step1_slaves = ('SEN_slaves_prep.py project.conf'); % slaves preparation
                        step2_slaves = ('SEN_splitting_slaves.py project.conf'); % slaves splitting & apply orbits
                        step3_slaves = ('SEN_coreg_ifg_topsar.py project.conf'); % coregistration & interferogram
                        step4_slaves = ('SEN_stamps_export.py project.conf'); % StaMPS export
                        step5_slaves = ('SEN_average_intensity.py project.conf'); % average instensity
                        step6_slaves = ('SEN_terrain_correction.py project.conf'); % terrain corrected coherence and lia

                            % CASES DEFENDING ON FIRST STEP
                            if isnumeric(first_step)
                                first_step_num = first_step;
                            elseif ischar(first_step)
                                first_step_num = str2double(first_step);
                            end
                            if first_step_num == 1
                                step_slaves_1 = [python space step1_slaves];
                                step_slaves_2 = [python space step2_slaves];
                                step_slaves_3 = [python space step3_slaves];
                                step_slaves_4 = [python space step4_slaves];
                                step_slaves_5 = [python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                end
                            elseif first_step_num == 2
                                step_slaves_1 = [dp python space step1_slaves];
                                step_slaves_2 = [python space step2_slaves];
                                step_slaves_3 = [python space step3_slaves];
                                step_slaves_4 = [python space step4_slaves];
                                step_slaves_5 = [python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                end
                            elseif first_step_num == 3
                                step_slaves_1 = [dp python space step1_slaves];
                                step_slaves_2 = [dp python space step2_slaves];
                                step_slaves_3 = [python space step3_slaves];
                                step_slaves_4 = [python space step4_slaves];
                                step_slaves_5 = [python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                end
                            elseif first_step_num == 4
                                step_slaves_1 = [dp python space step1_slaves];
                                step_slaves_2 = [dp python space step2_slaves];
                                step_slaves_3 = [dp python space step3_slaves];
                                step_slaves_4 = [python space step4_slaves];
                                step_slaves_5 = [python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                end
                             elseif first_step_num == 5
                                step_slaves_1 = [dp python space step1_slaves];
                                step_slaves_2 = [dp python space step2_slaves];
                                step_slaves_3 = [dp python space step3_slaves];
                                step_slaves_4 = [dp python space step4_slaves];
                                step_slaves_5 = [python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                end
                            elseif first_step_num == 6
                                step_slaves_1 = [dp python space step1_slaves];
                                step_slaves_2 = [dp python space step2_slaves];
                                step_slaves_3 = [dp python space step3_slaves];
                                step_slaves_4 = [dp python space step4_slaves];
                                step_slaves_5 = [dp python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
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
                            fprintf(f_snap2stamps_slaves,'#!/bin/bash \n');
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
                            fprintf(f_snap2stamps_slaves,'#!/bin/bash \n');
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
                            system(chmod);
                            phase_preprocessing_beta.runCommandLive(app, path_2_slaves, ...
                                'Slave processing pipeline', 18, 92, first_step_num, 6);
                        elseif ismac
                            path_2_slaves = strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh');
                            chmod_in = 'chmod +x';
                            chmod = [chmod_in space path_2_slaves];
                            % Executed without opening Terminal by runCommandLive.
                            system(chmod);
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
                                fprintf(f_average_intensity,'%s\r\n',[python space step5_slaves]);
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
                                fprintf(f_average_intensity,'#!/bin/bash \n');
                                fprintf(f_average_intensity,'cd "%s"\n',path_1_slaves);
                                fprintf(f_average_intensity,'%s\n',[python space step5_slaves]);
                                fprintf(f_average_intensity,'sleep 5 \n');
                                fprintf(f_average_intensity,'cd "%s"\n',path_1_download);
                                fprintf(f_average_intensity,'touch dummyAverageIntensity.txt \n');
                                f_average_intensity = fclose(f_average_intensity);
                                path_2_average_intensity = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_update_average_intensity.sh'));
                                chmod_in = ('chmod +x');
                                chmod = [chmod_in space path_2_average_intensity];
                                xterm = ('sudo xterm -e');
                                system(chmod);
                                phase_preprocessing_beta.runCommandLive(app, path_2_average_intensity, ...
                                    'Full-stack average intensity', 78, 92, 5, 5);
                            elseif ismac
                                f_average_intensity = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_update_average_intensity.sh'),'w');
                                j = 0;
                                while f_average_intensity == -1 && j < 20
                                    f_average_intensity = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_update_average_intensity.sh'),'w');
                                    j = j+1;
                                end
                                fprintf(f_average_intensity,'#!/bin/bash \n');
                                fprintf(f_average_intensity,'cd "%s"\n',path_1_slaves);
                                fprintf(f_average_intensity,'%s\n',[python space step5_slaves]);
                                fprintf(f_average_intensity,'sleep 5 \n');
                                fprintf(f_average_intensity,'cd "%s"\n',path_1_download);
                                fprintf(f_average_intensity,'touch dummyAverageIntensity.txt \n');
                                f_average_intensity = fclose(f_average_intensity);
                                path_2_average_intensity = strcat(project_path_full, '/snap2stamps/bin/snap2stamps_update_average_intensity.sh');
                                chmod_in = 'chmod +x';
                                chmod = [chmod_in space path_2_average_intensity];
                                % Executed without opening Terminal by runCommandLive.
                                system(chmod);
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

                        project_parent_path_full = prep_folder;
                        stamps_folder_full = fullfile(project_parent_path_full, stamps_folder);
                        if ~isfolder(stamps_folder_full)
                            mkdir(stamps_folder_full);
                            updateOutput(app, ['Created StaMPS dataset folder: ' stamps_folder_full]);
                        end
                        % The beta launches the editable StaMPS module from its
                        % canonical installation. No PHASE_StaMPS.mlapp is copied
                        % or required at runtime.
                        stamps_diagnostic = fullfile(project_path_full, 'diagnose_PHASE_StaMPS.m');
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

                        % OPEN PHASE_StaMPS_beta (cross-platform). We change MATLAB
                        % current directory to the StaMPS folder so the relative
                        % paths inside PHASE_StaMPS (./input_StaMPS.mat, INSAR_*/,
                        % diff0/, ...) resolve correctly, then launch the app.
                        stamps_app_full = fullfile(project_parent_path_full, stamps_folder);
                        stamps_app_file = fullfile(project_path_full, 'PHASE_StaMPS_beta.m');
                        updateOutput(app, ['Preprocessing completed. StaMPS dataset folder: ' stamps_app_full]);
                        choice = 'Open now';
                        if isfile(dst_input_mat)
                            updateOutput(app, ['Opening PHASE_StaMPS_beta with the existing configuration in: ' stamps_app_full]);
                        else
                            updateOutput(app, ['Opening PHASE_StaMPS_beta to create the initial configuration in: ' stamps_app_full]);
                        end
                        if strcmp(choice, 'Open now')
                            if isfile(stamps_app_file)
                                updateOutput(app, ['Launching the canonical PHASE_StaMPS_beta runtime for: ' stamps_app_full]);
                                try
                                    phase_preprocessing_beta.launchStampsBeta(stamps_app_file, stamps_app_full);
                                catch ex
                                    updateOutput(app, ['Failed to launch PHASE_StaMPS: ' ex.message]);
                                    updateOutput(app, ['Open it manually with: PHASE_StaMPS_beta(''' stamps_app_full ''')']);
                                end
                            else
                                updateOutput(app, ['PHASE_StaMPS_beta.m not found in ' project_path_full ', please open it manually.']);
                            end
                        else
                            updateOutput(app, ['You can open PHASE_StaMPS later from: ' stamps_app_full]);
                        end

                        %% ----------------------------------------------------

                        % Display a message when the script is done
                        updateOutput(app, 'Script SEN_Preprocessing.m completed.');

                        % Change lamp color to green (success)
                        app.PreprocessingstatusLamp.Color = [0, 1, 0]; % Green

                    catch ME

                        % If an error occurs
                        if strcmp(ME.identifier,'PHASE:ProcessingStopped')
                            updateOutput(app, 'Preprocessing was force-stopped by the user.');
                        else
                            updateOutput(app, 'An error occurred during script execution. Please check each step log file!');
                        end

                        % Change lamp color to red (error)
                        app.PreprocessingstatusLamp.Color = [1, 0, 0]; % Red

                        % Show the actual error
                        rethrow(ME)

                    end

                    condition = 1;
                end

            elseif app.constellation == "CSK"

                condition = 0;

                while condition == 0

                    app.MessagesTextArea.Value = '';

                    % Display a message in the UIAxes
                    updateOutput(app, 'Button pressed! Running the script...');
                    drawnow;

                    % Change lamp color to yellow (running)
                    app.PreprocessingstatusLamp_2.Color = [1, 1, 0]; % Yellow

                    try

                        %% COSMO-SKYMED PROCESSING ----------------------------

                        pwd;
                        prep_folder = pwd;
                        addpath(genpath(prep_folder));

                        %% ------------------ IMPORT OF THE REQUIRED VARIABLES --------------------

                        % READ OF THE INPUT VARIABLES
                        par = filesep;

                        load(strcat('.', par, 'PHASE_Preprocessing', par, 'input_preprocessing.mat'), 'python', ...
                                                'master_date', 'auto_master', 'master_processing', 'lon_min', ...
                                                'lat_min', 'lon_max', 'lat_max', 'slaves_removal', 'dem_name', 'dem_file', ...
                                                'first_step', 'num_gcp', 'coherence_tc', 'epsg_code', 'gptbin_path', 'cpu', 'cache');

                        % PROJECT FOLDER
                        project_path_full = strcat(prep_folder, par, 'PHASE_Preprocessing');
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
                            graphsfolder = strcat('GRAPHSFOLDER=', project_path_full, par, 'snap2stamps', par, 'graphs');

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
                            step_selector = ('CSK_master_selector.py project_master.conf');
                            step_selector_cmd = [python space step_selector];
                            step1_master = ('CSK_subset_master.py project_master.conf');
                            step_master_2 = [python space step1_master];

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
                                fprintf(f_snap2stamps_master,'#!/bin/bash \n');
                                fprintf(f_snap2stamps_master,'cd "%s"\n',path_1_master);
                                fprintf(f_snap2stamps_master,'%s\n',step_selector_cmd);
                                fprintf(f_snap2stamps_master,'%s\n',step_master_2);
                                fclose(f_snap2stamps_master);

                                path_2_master = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'));
                                chmod = ['chmod +x' space path_2_master];
                                system(chmod);
                                phase_preprocessing_beta.runCommandLive(app, path_2_master, ...
                                    'Master selection and preparation', 3, 18, 1, 1);

                            elseif ismac
                                f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                j = 0;
                                while f_snap2stamps_master == -1 && j < 20
                                    f_snap2stamps_master = fopen(strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'),'w');
                                    j = j+1;
                                end
                                fprintf(f_snap2stamps_master,'#!/bin/bash \n');
                                fprintf(f_snap2stamps_master,'cd "%s"\n',path_1_master);
                                fprintf(f_snap2stamps_master,'%s\n',step_selector_cmd);
                                fprintf(f_snap2stamps_master,'%s\n',step_master_2);
                                fclose(f_snap2stamps_master);

                                path_2_master = (strcat(project_path_full, '/snap2stamps/bin/snap2stamps_master.sh'));
                                chmod = ['chmod +x' space path_2_master];
                                system(chmod);
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
                        graphsfolder = strcat('GRAPHSFOLDER=', project_path_full, par, 'snap2stamps', par, 'graphs'); % path of the graphs folder of snap2stamps
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
                        dp = ('::');
                        step1_slaves = ('CSK_slaves_prep.py project.conf'); % slaves preparation
                        step2_slaves = ('CSK_subset_slaves.py project.conf'); % slaves splitting & apply orbits
                        step3_slaves = ('CSK_coreg_ifg.py project.conf'); % coregistration & interferogram
                        step4_slaves = ('CSK_stamps_export.py project.conf'); % StaMPS export
                        step5_slaves = ('CSK_average_intensity.py project.conf'); % terrain corrected coherence and lia
                        step6_slaves = ('CSK_terrain_correction.py project.conf'); % terrain corrected coherence and lia

                            % CASES DEFENDING ON FIRST STEP
                            if isnumeric(first_step)
                                first_step_num = first_step;
                            elseif ischar(first_step)
                                first_step_num = str2double(first_step);
                            end
                            if first_step_num == 1
                                step_slaves_1 = [python space step1_slaves];
                                step_slaves_2 = [python space step2_slaves];
                                step_slaves_3 = [python space step3_slaves];
                                step_slaves_4 = [python space step4_slaves];
                                step_slaves_5 = [python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                end
                            elseif first_step_num == 2
                                step_slaves_1 = [dp python space step1_slaves];
                                step_slaves_2 = [python space step2_slaves];
                                step_slaves_3 = [python space step3_slaves];
                                step_slaves_4 = [python space step4_slaves];
                                step_slaves_5 = [python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                end
                            elseif first_step_num == 3
                                step_slaves_1 = [dp python space step1_slaves];
                                step_slaves_2 = [dp python space step2_slaves];
                                step_slaves_3 = [python space step3_slaves];
                                step_slaves_4 = [python space step4_slaves];
                                step_slaves_5 = [python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                end
                            elseif first_step_num == 4
                                step_slaves_1 = [dp python space step1_slaves];
                                step_slaves_2 = [dp python space step2_slaves];
                                step_slaves_3 = [dp python space step3_slaves];
                                step_slaves_4 = [python space step4_slaves];
                                step_slaves_5 = [python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                end
                            elseif first_step_num == 5
                                step_slaves_1 = [dp python space step1_slaves];
                                step_slaves_2 = [dp python space step2_slaves];
                                step_slaves_3 = [dp python space step3_slaves];
                                step_slaves_4 = [dp python space step4_slaves];
                                step_slaves_5 = [python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                end
                            elseif first_step_num == 6
                                step_slaves_1 = [dp python space step1_slaves];
                                step_slaves_2 = [dp python space step2_slaves];
                                step_slaves_3 = [dp python space step3_slaves];
                                step_slaves_4 = [dp python space step4_slaves];
                                step_slaves_5 = [dp python space step5_slaves];
                                if coherence_tc == 0
                                    step_slaves_6 = [python space step6_slaves];
                                else
                                    updateOutput(app, ['You are running just the terrain correction for coherence and LIA bands, but you have set tc_coherence = 1. ' ...
                                        'Please change it to 0 to perform this step.'])
                                end
                            end

                        if update_processed_run
                            step_slaves_5 = [dp python space step5_slaves];
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
                            fprintf(f_snap2stamps_slaves,'#!/bin/bash \n');
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
                            fprintf(f_snap2stamps_slaves,'#!/bin/bash \n');
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
                            system(chmod);
                            phase_preprocessing_beta.runCommandLive(app, path_2_slaves, ...
                                'Slave processing pipeline', 18, 92, first_step_num, 6);
                        elseif ismac
                            path_2_slaves = strcat(project_path_full, '/snap2stamps/bin/snap2stamps_slaves.sh');
                            chmod_in = 'chmod +x';
                            chmod = [chmod_in space path_2_slaves];
                            % Executed without opening Terminal by runCommandLive.
                            system(chmod);
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

                        project_parent_path_full = prep_folder;
                        stamps_folder_full = fullfile(project_parent_path_full, stamps_folder);
                        if ~isfolder(stamps_folder_full)
                            mkdir(stamps_folder_full);
                            updateOutput(app, ['Created StaMPS dataset folder: ' stamps_folder_full]);
                        end
                        % The beta launches the editable StaMPS module from its
                        % canonical installation. No PHASE_StaMPS.mlapp is copied
                        % or required at runtime.
                        stamps_diagnostic = fullfile(project_path_full, 'diagnose_PHASE_StaMPS.m');
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

                        % OPEN PHASE_StaMPS_beta (cross-platform). We change MATLAB
                        % current directory to the StaMPS folder so the relative
                        % paths inside PHASE_StaMPS (./input_StaMPS.mat, INSAR_*/,
                        % diff0/, ...) resolve correctly, then launch the app.
                        stamps_app_full = fullfile(project_parent_path_full, stamps_folder);
                        stamps_app_file = fullfile(project_path_full, 'PHASE_StaMPS_beta.m');
                        updateOutput(app, ['Preprocessing completed. StaMPS dataset folder: ' stamps_app_full]);
                        choice = 'Open now';
                        if isfile(dst_input_mat)
                            updateOutput(app, ['Opening PHASE_StaMPS_beta with the existing configuration in: ' stamps_app_full]);
                        else
                            updateOutput(app, ['Opening PHASE_StaMPS_beta to create the initial configuration in: ' stamps_app_full]);
                        end
                        if strcmp(choice, 'Open now')
                            if isfile(stamps_app_file)
                                updateOutput(app, ['Launching the canonical PHASE_StaMPS_beta runtime for: ' stamps_app_full]);
                                try
                                    phase_preprocessing_beta.launchStampsBeta(stamps_app_file, stamps_app_full);
                                catch ex
                                    updateOutput(app, ['Failed to launch PHASE_StaMPS: ' ex.message]);
                                    updateOutput(app, ['Open it manually with: PHASE_StaMPS_beta(''' stamps_app_full ''')']);
                                end
                            else
                                updateOutput(app, ['PHASE_StaMPS_beta.m not found in ' project_path_full ', please open it manually.']);
                            end
                        else
                            updateOutput(app, ['You can open PHASE_StaMPS later from: ' stamps_app_full]);
                        end

                        %% ----------------------------------------------------

                        % Display a message when the script is done
                        updateOutput(app, 'Script CSK_Preprocessing.m completed.');

                        % Change lamp color to green (success)
                        app.PreprocessingstatusLamp_2.Color = [0, 1, 0]; % Green

                    catch ME

                        % If an error occurs
                        if strcmp(ME.identifier,'PHASE:ProcessingStopped')
                            updateOutput(app, 'Preprocessing was force-stopped by the user.');
                        else
                            updateOutput(app, 'An error occurred during script execution. Please check each step log file!');
                        end

                        % Change lamp color to red (error)
                        app.PreprocessingstatusLamp_2.Color = [1, 0, 0]; % Red

                        % Show the actual error
                        rethrow(ME)

                    end

                    condition = 1;

                end

            end

        end

        % Button pushed function: LoadButton, LoadButton_2
        function LoadButtonPushed(app, event)

            filename = './PHASE_Preprocessing/input_preprocessing.mat';  % Specify the filename

            if exist(filename, 'file') == 2
                data = load(filename);
                app.constellation = data.constellation;

                if app.constellation == "SEN"
                    % Update SEN-related properties
                    if data.python == "python"   % python
                        app.PythonEnvironmentDropDown.Value = 'python';
                    elseif data.python == "python2"
                        app.PythonEnvironmentDropDown.Value = 'python2';
                    elseif data.python == "python2.7"
                        app.PythonEnvironmentDropDown.Value = 'python2.7';
                    else
                        app.PythonEnvironmentDropDown.Value = 'Other';
                        app.PythonEnvironmentLabel.Visible = 'on';
                        app.CustomPythonEnvironmentEditField.Visible = 'on';
                    end
                    app.CustomPythonEnvironmentEditField.Value = data.python;
                    date_txt = data.master_date;   % master date
                    date_dt = datetime(date_txt, 'InputFormat', 'yyyyMMdd', 'Format', 'dd-MMM-uuuu');
                    app.MasterdateDatePicker.Value = date_dt;
                    if isfield(data, 'auto_master')
                        app.AutoMasterCheckBox.Value = logical(data.auto_master);
                        app.auto_master_SEN = data.auto_master;
                        if data.auto_master == 1
                            app.MasterdateDatePicker.Enable = 'off';
                        else
                            app.MasterdateDatePicker.Enable = 'on';
                        end
                    end
                    if data.master_processing == 0   % master processing
                        app.MasterprocessingCheckBox.Value = true;
                    elseif data.master_processing == 1
                        app.MasterprocessingCheckBox.Value = false;
                    end
                    app.PolarisationDropDown.Value = data.polarisation;
                    app.MinlongitudeEditField.Value = data.lon_min;
                    app.MinlatitudeEditField.Value = data.lat_min;
                    app.MaxlongitudeEditField.Value = data.lon_max;
                    app.MaxlatitudeEditField.Value = data.lat_max;
                    % Move the red box on the Sentinel-1 Map
                    if isvalid(app.roi_SEN)
                        w = data.lon_max - data.lon_min;
                        h = data.lat_max - data.lat_min;
                        % GeoAxes expects [LatMin, LonMin, LatHeight, LonWidth]
                        app.roi_SEN.Position = [data.lat_min, data.lon_min, h, w];
                    end
                    if data.slaves_removal == 0   % slaves removal
                        app.SlavesremovalafterprocessingCheckBox.Value = true;
                    elseif data.slaves_removal == 1
                        app.SlavesremovalafterprocessingCheckBox.Value = false;
                    end
                    app.DEMinterferogramDropDown.Value = data.dem_name;
                    app.DEMifgpathEditField.Value = data.dem_file;
                    app.DEMcoregistrationDropDown.Value = data.dem_name_coreg;
                    app.DEMcoregpathEditField.Value = data.dem_file_coreg;
                    app.DEMresamplingmethodDropDown.Value = data.dem_resampling;
                    app.FirststepDropDown.Value = num2str(data.first_step);
                    if data.coherence_tc == 0   % coherence tc
                        app.TerraincorrectedCoherenceandLIACheckBox.Value = true;
                    elseif data.coherence_tc == 1
                        app.TerraincorrectedCoherenceandLIACheckBox.Value = false;
                    end
                    app.EPSGcodeEditField.Value = data.epsg_code;
                    app.PathEditField.Value = data.gptbin_path;
                    app.CPUEditField.Value = data.cpu;
                    app.CacheEditField.Value = data.cache;
                    drawnow;

                    % Set the lamp color to green to indicate successful loading
                    app.StatusLamp_3.Color = [0, 1, 0];  % Green

                    % Update also the variables for the saving
                    app.python_SEN = app.CustomPythonEnvironmentEditField.Value;
                    app.master_date_SEN = data.master_date;
                    app.auto_master_SEN = data.auto_master;
                    app.master_processing_SEN = data.master_processing;
                    app.polarisation_SEN = app.PolarisationDropDown.Value;
                    app.lon_min_SEN = app.MinlongitudeEditField.Value;
                    app.lat_min_SEN = app.MinlatitudeEditField.Value;
                    app.lon_max_SEN = app.MaxlongitudeEditField.Value;
                    app.lat_max_SEN = app.MaxlatitudeEditField.Value;
                    app.slaves_removal_SEN = data.slaves_removal;
                    app.dem_name_SEN = app.DEMinterferogramDropDown.Value;
                    app.dem_file_SEN = app.DEMifgpathEditField.Value;
                    app.dem_name_coreg_SEN = app.DEMcoregistrationDropDown.Value;
                    app.dem_file_coreg_SEN = app.DEMcoregpathEditField.Value;
                    app.dem_resampling_SEN = app.DEMresamplingmethodDropDown.Value;
                    app.first_step_SEN = app.FirststepDropDown.Value;
                    app.coherence_tc_SEN = data.coherence_tc;
                    app.epsg_code_SEN = app.EPSGcodeEditField.Value;
                    app.gptbin_path_SEN = app.PathEditField.Value;
                    app.cpu_SEN = app.CPUEditField.Value;
                    app.cache_SEN = app.CacheEditField.Value;

                elseif app.constellation == "CSK"
                    % Update CSK-related properties
                    if data.python == "python"   % python
                        app.PythonEnvironmentDropDown_2.Value = 'python';
                    elseif data.python == "python2"
                        app.PythonEnvironmentDropDown_2.Value = 'python2';
                    elseif data.python == "python2.7"
                        app.PythonEnvironmentDropDown_2.Value = 'python2.7';
                    else
                        app.PythonEnvironmentDropDown_2.Value = 'Other';
                        app.CustomPythonEnvironmentLabel.Visible = 'off';
                        app.CustomPythonEnvironmentEditField_2.Visible = 'off';
                    end
                    app.CustomPythonEnvironmentEditField_2.Value = data.python;
                    date_txt = data.master_date;   % master date
                    date_dt = datetime(date_txt, 'InputFormat', 'yyyyMMdd', 'Format', 'dd-MMM-uuuu');
                    app.MasterdateDatePicker_2.Value = date_dt;
                    if isfield(data, 'auto_master')
                        app.AutoMasterCheckBox_2.Value = logical(data.auto_master);
                        app.auto_master_CSK = data.auto_master;
                        if data.auto_master == 1
                            app.MasterdateDatePicker_2.Enable = 'off';
                        else
                            app.MasterdateDatePicker_2.Enable = 'on';
                        end
                    end
                    if data.master_processing == 0   % master processing
                        app.MasterprocessingCheckBox_2.Value = true;
                    elseif data.master_processing == 1
                        app.MasterprocessingCheckBox_2.Value = false;
                    end
                    app.MinlongitudeEditField_2.Value = data.lon_min;
                    app.MinlatitudeEditField_2.Value = data.lat_min;
                    app.MaxlongitudeEditField_2.Value = data.lon_max;
                    app.MaxlatitudeEditField_2.Value = data.lat_max;
                    % Move the red box on the Cosmo-SkyMed Map
                    if isvalid(app.roi_CSK)
                        w = data.lon_max - data.lon_min;
                        h = data.lat_max - data.lat_min;
                        % GeoAxes expects [LatMin, LonMin, LatHeight, LonWidth]
                        app.roi_CSK.Position = [data.lat_min, data.lon_min, h, w];
                    end
                    if data.slaves_removal == 0   % slaves removal
                        app.SlavesremovalafterprocessingCheckBox_2.Value = true;
                    elseif data.slaves_removal == 1
                        app.SlavesremovalafterprocessingCheckBox_2.Value = false;
                    end
                    app.DEMinterferogramDropDown_2.Value = data.dem_name;
                    app.DEMifgpathEditField_2.Value = data.dem_file;
                    app.FirststepDropDown_2.Value = num2str(data.first_step);
                    app.CoregistrationGCPsnumberEditField.Value = data.num_gcp;
                    if data.coherence_tc == 0   % coherence tc
                        app.TerraincorrectedCoherenceandLIACheckBox_2.Value = true;
                    elseif data.coherence_tc == 1
                        app.TerraincorrectedCoherenceandLIACheckBox_2.Value = false;
                    end
                    app.EPSGcodeEditField_2.Value = data.epsg_code;
                    app.PathEditField_2.Value = data.gptbin_path;
                    app.CPUEditField_2.Value = data.cpu;
                    app.CacheEditField_2.Value = data.cache;
                    drawnow;

                    % Set the lamp color to green to indicate successful loading
                    app.StatusLamp_4.Color = [0, 1, 0];  % Green

                    % Update also the variables for the saving
                    app.python_CSK = app.CustomPythonEnvironmentEditField_2.Value;
                    app.master_date_CSK = data.master_date;
                    app.auto_master_CSK = data.auto_master;
                    app.master_processing_CSK = data.master_processing;
                    app.lon_min_CSK = app.MinlongitudeEditField_2.Value;
                    app.lat_min_CSK = app.MinlatitudeEditField_2.Value;
                    app.lon_max_CSK = app.MaxlongitudeEditField_2.Value;
                    app.lat_max_CSK = app.MaxlatitudeEditField_2.Value;
                    app.slaves_removal_CSK = data.slaves_removal;
                    app.dem_name_CSK = app.DEMinterferogramDropDown_2.Value;
                    app.dem_file_CSK = app.DEMifgpathEditField_2.Value;
                    app.num_gcp_CSK = app.CoregistrationGCPsnumberEditField.Value;
                    app.first_step_CSK = app.FirststepDropDown_2.Value;
                    app.coherence_tc_CSK = data.coherence_tc;
                    app.epsg_code_CSK = app.EPSGcodeEditField_2.Value;
                    app.gptbin_path_CSK = app.PathEditField_2.Value;
                    app.cpu_CSK = app.CPUEditField_2.Value;
                    app.cache_CSK = app.CacheEditField_2.Value;

                end


            else
                % Handle the case when the file doesn't exist
                % Set the lamp color to red to indicate an error
                app.StatusLamp_3.Color = [1, 0, 0];  % Red
                app.StatusLamp_4.Color = [1, 0, 0];  % Red
            end

        end

        % Selection change function: TabGroup
        function TabGroupSelectionChanged(app, event)

            % Reset lamp colors to red
            app.StatusLamp.Color = [1, 0, 0]; % Red color for Save lamp
            app.StatusLamp_3.Color = [1, 0, 0]; % Red color for Load lamp

            if isequal(app.TabGroup.SelectedTab, app.ImagesTab_SEN)
                app.ImportedSENImagesTable.Data = cell(0, 3);
                drawnow;
                app.refreshImportedSENImagesTable();
            end

            % % Controlla se il tab selezionato è quello dell'AOI
            % if event.NewValue == app.AOITab_2
            %     app.initializeSENMap();
            % end

        end

        % Selection change function: TabGroup2
        function TabGroup2SelectionChanged(app, event)

            % Reset lamp colors to red
            app.StatusLamp_2.Color = [1, 0, 0]; % Red color for Save lamp
            app.StatusLamp_4.Color = [1, 0, 0]; % Red color for Load lamp

            % % Controlla se il tab selezionato è quello dell'AOI
            % if event.NewValue == app.AOITab
            %     app.initializeCSKMap();
            % end

        end

        % Button pushed function: BrowseButton
        function BrowseButtonPushed(app, event)

            % Open a file dialog for the user to select the external DEM file
            [filename, filepath] = uigetfile('*.tif', 'Select .tif DEM file');

            % Check if the user selected a file (clicked "OK" in the file dialog)
            if isequal(filename, 0)
                % User canceled the file selection, do nothing
            else
                % User selected a file, update the EditField with the selected file path
                selectedFilePath = fullfile(filepath, filename);
                app.DEMifgpathEditField.Value = selectedFilePath;
                app.dem_file_SEN = selectedFilePath;
            end

        end

        % Button pushed function: BrowseButton_2
        function BrowseButton_2Pushed(app, event)

            % Open a file dialog for the user to select the external DEM file
            [filename, filepath] = uigetfile('*.tif', 'Select .tif DEM file');

            % Check if the user selected a file (clicked "OK" in the file dialog)
            if isequal(filename, 0)
                % User canceled the file selection, do nothing
            else
                % User selected a file, update the EditField with the selected file path
                selectedFilePath = fullfile(filepath, filename);
                app.DEMcoregpathEditField.Value = selectedFilePath;
                app.dem_file_coreg_SEN = selectedFilePath;
            end

        end

        % Button pushed function: BrowseButton_3
        function BrowseButton_3Pushed(app, event)

            % Open a file dialog for the user to select the external DEM file
            [filename, filepath] = uigetfile('*.tif', 'Select .tif DEM file');

            % Check if the user selected a file (clicked "OK" in the file dialog)
            if isequal(filename, 0)
                % User canceled the file selection, do nothing
            else
                % User selected a file, update the EditField with the selected file path
                selectedFilePath = fullfile(filepath, filename);
                app.DEMifgpathEditField_2.Value = selectedFilePath;
                app.dem_file_CSK = selectedFilePath;
            end

        end

        % Value changed function: PythonEnvironmentDropDown
        function PythonEnvironmentDropDownValueChanged(app, event)
            value = app.PythonEnvironmentDropDown.Value;

            % Check if "Other" is selected
            if strcmp(value, 'Other')
                app.PythonEnvironmentLabel.Visible = 'on'; % Show the Label
                app.CustomPythonEnvironmentEditField.Visible = 'on'; % Show the EditField
                app.python_SEN = app.CustomPythonEnvironmentEditField.Value;
            else
                app.PythonEnvironmentLabel.Visible = 'off'; % Hide the Label
                app.CustomPythonEnvironmentEditField.Visible = 'off'; % Hide the EditField
                app.python_SEN = value;
            end

        end

        % Value changed function: PythonEnvironmentDropDown_2
        function PythonEnvironmentDropDown_2ValueChanged(app, event)
            value = app.PythonEnvironmentDropDown_2.Value;

            % Check if "Other" is selected
            if strcmp(value, 'Other')
                app.CustomPythonEnvironmentLabel.Visible = 'on'; % Show the Label
                app.CustomPythonEnvironmentEditField_2.Visible = 'on'; % Show the EditField
                app.python_CSK = app.CustomPythonEnvironmentEditField_2.Value;
            else
                app.CustomPythonEnvironmentLabel.Visible = 'off'; % Hide the Label
                app.CustomPythonEnvironmentEditField_2.Visible = 'off'; % Hide the EditField
                app.python_CSK = value;
            end

        end

        % Button pushed function: StopButton, StopButton_2
        function StopButtonPushed(app, event)

            % Set a flag to indicate that the button has been pressed
            app.StopFlag = true;

        end

        % Value changed function: AutoMasterCheckBox
        function AutoMasterCheckBoxValueChanged(app, event)
            value = app.AutoMasterCheckBox.Value;
            if value
                app.auto_master_SEN = 1;
                app.MasterdateDatePicker.Enable = 'off';
            else
                app.auto_master_SEN = 0;
                app.MasterdateDatePicker.Enable = 'on';
            end
        end

        % Value changed function: AutoMasterCheckBox_2
        function AutoMasterCheckBox_2ValueChanged(app, event)
            value = app.AutoMasterCheckBox_2.Value;
            if value
                app.auto_master_CSK = 1;
                app.MasterdateDatePicker_2.Enable = 'off';
            else
                app.auto_master_CSK = 0;
                app.MasterdateDatePicker_2.Enable = 'on';
            end
        end

        % Button pushed function: DrawCSKAOIButton
        function DrawCSKAOIButtonPushed(app, event)
            app.initializeCSKMap();
            updateOutput(app, 'Map and AOI reset to default footprint.');
        end

        % Button pushed function: ResetSENMapButton
        function ResetSENMapButtonPushed(app, event)
            app.initializeSENMap();
            updateOutput(app, 'Sentinel-1 Map and AOI reset to default footprint.');
        end

        % Callback function
        function SaveParametersButtonPushed(app, event)


            params = getDownloaderParams(app);

            outputFolder = fullfile(pwd, "downloadasf");

            if ~exist(outputFolder, "dir")
                mkdir(outputFolder);
            end

            outputFile = fullfile(outputFolder, "search_request.json");

            jsonText = jsonencode(params, PrettyPrint=true);

            fid = fopen(outputFile, "w");



            fprintf(fid, "%s", jsonText);
            fclose(fid);

            disp("Parameters saved to download/search_request.json");

            disp(jsonText);
        end

        % Button pushed function: ClearPathFrameButton
        function ClearPathFrameButtonPushed(app, event)
            app.PathStartEditField.Value = 0;
            app.PathEndEditField.Value = 0;
            app.FrameStartEditField.Value = 0;
            app.FrameEndEditField.Value = 0;
        end

        % Button pushed function: ShowFiltersButton
        function ShowFiltersButtonPushed(app, event)
            if strcmp(app.FilterPanel.Visible, "on")
                app.FilterPanel.Visible = "off";
                app.ShowFiltersButton.Text = "Show Filters";
            else
                app.FilterPanel.Visible = "on";
                app.ShowFiltersButton.Text = "Hide Filters";
            end
        end

        % Button pushed function: ClearAOIButton
        function ClearAOIButtonPushed(app, event)
            clearDownloaderAOIShape(app);

            app.DownloaderAOIType = "";
            app.DownloaderAOICorners = [];
            app.DownloaderPolygonCoords = [];

            app.MinLongitudeEditField.Value = -180;
            app.MaxLongitudeEditField.Value = 180;
            app.MinLatitudeEditField.Value = -90;
            app.MaxLatitudeEditField.Value = 90;
        end

        % Button pushed function: DrawRectangleButton
        function DrawRectangleButtonPushed(app, event)
            ClearAOIButtonPushed(app);

            app.roi_Downloader = drawrectangle(app.DownloaderGeoAxes, ...
                'Color', 'r', ...
                'FaceAlpha', 0.2);

            app.updateDownloaderRectangleCoords(app.roi_Downloader);

            addlistener(app.roi_Downloader, ...
                'ROIMoved', ...
                @(src, event) app.updateDownloaderRectangleCoords(src));

        end

        % Button pushed function: SearchASFButton
        function SearchASFButtonPushed(app, event)

            if ~hasDownloaderAOI(app)
                uialert(app.UIFigure, ...
                    "Please draw or select an Area of Interest before searching.", ...
                    "Missing AOI", ...
                    "Icon", "warning");
                return;
            end
            SaveParametersButtonPushed(app, event);
            disp("Fresh search_request.json:");
            disp(fileread(fullfile(pwd, "downloadasf", "search_request.json")));
            disp("Search ASF button pressed");

            appPath = phase_preprocessing_beta.projectRoot();
            controllerPath = fullfile(appPath, "downloadasf", "controller.py");

            pythonEnv = pyenv;
            pythonExe = string(pythonEnv.Executable);

            command = """" + pythonExe + """ """ + controllerPath + """ search";

            disp("Command being run:");
            disp(command);

            [status, cmdout] = system(command);

            disp("Python status:");
            disp(status);

            disp("Python output:");
            disp(cmdout);

            if status ~= 0
                uialert(app.UIFigure, cmdout, "Search Failed", "Icon", "error");
                return;
            end
            backupSearchDownloadFiles(app);
            loadSearchPreview(app);
            previewRequested = showSearchPreviewChoice(app);

            if previewRequested
                app.PreviewPanel.Visible = "on";
                disp("Preview requested.");
            else
                disp("Preview skipped.");
            end
        end

        % Button pushed function: ResetallfiltersButton
        function ResetallfiltersButtonPushed(app, event)
            initializeParams(app);
        end

        % Button pushed function: DrawPolygonButton
        function DrawPolygonButtonPushed(app, event)
            ClearAOIButtonPushed(app);

            app.roi_Downloader = drawpolygon(app.DownloaderGeoAxes, ...
                'Color', 'r', ...
                'FaceAlpha', 0.2);

            updateDownloaderPolygonCoords(app, app.roi_Downloader);

            addlistener(app.roi_Downloader, ...
                'ROIMoved', ...
                @(src, event) updateDownloaderPolygonCoords(app, src));
        end

        % Value changed function: MaxLatitudeEditField
        function MaxLatitudeEditFieldValueChanged(app, event)
            value = app.MaxLatitudeEditField.Value;
            AOIBoundaryValueChanged(app, event);
        end

        % Value changed function: MinLatitudeEditField
        function MinLatitudeEditFieldValueChanged(app, event)
            value = app.MinLatitudeEditField.Value;
            AOIBoundaryValueChanged(app, event);
        end

        % Value changed function: MaxLongitudeEditField
        function MaxLongitudeEditFieldValueChanged(app, event)
            value = app.MaxLongitudeEditField.Value;
            AOIBoundaryValueChanged(app, event);
        end

        % Value changed function: MinLongitudeEditField
        function MinLongitudeEditFieldValueChanged(app, event)
            value = app.MinLongitudeEditField.Value;
            AOIBoundaryValueChanged(app, event);
        end

        % Button pushed function: LoginButton
        function LoginButtonPushed(app, event)
            username = strtrim(string(app.EarthdataUsernameEditField.Value));
            password = strtrim(string(app.EarthdataPasswordEditField.Value));

            app.LoginFeedbackLabel.Text = "Logging in...";

            if username == "" || password == ""
                app.LoginFeedbackLabel.Text = ...
                    "Username and password are required.";
                return;
            end

            loginData.username = username;
            loginData.password = password;

            outputFolder = fullfile(pwd, "downloadasf");

            if ~exist(outputFolder, "dir")
                mkdir(outputFolder);
            end

            outputFile = fullfile(outputFolder, "login_request.json");

            fid = fopen(outputFile, "w");
            fprintf(fid, "%s", jsonencode(loginData, PrettyPrint=true));
            fclose(fid);

            app.LoginFeedbackLabel.Text = "Checking login...";

            % backedn logic
            appPath = phase_preprocessing_beta.projectRoot();
            controllerPath = fullfile(appPath, "downloadasf", "controller.py");

            pythonEnv = pyenv;
            pythonExe = string(pythonEnv.Executable);

            command = """" + pythonExe + """ """ + controllerPath + """ login";

            disp("Running login command:");
            disp(command);

            [status, cmdout] = system(command);

            disp("Login backend status:");
            disp(status);
            disp("Login backend output:");
            disp(cmdout);

            if status ~= 0
                app.LoginFeedbackLabel.Text = "Login backend failed.";
                return;
            end

            resultFile = fullfile(appPath, "downloadasf", "login_result.json");

            if ~exist(resultFile, "file")
                app.LoginFeedbackLabel.Text = "Missing login_result.json.";
                return;
            end

            result = jsondecode(fileread(resultFile));

            loginSuccess = strcmp(string(result.status), "success");

            if loginSuccess

                app.LoginFeedbackLabel.Text = "";

                app.EarthdataUsernameEditField.Visible = "off";
                app.EarthdataPasswordEditField.Visible = "off";
                app.LoginButton.Visible = "off";

                app.SignedInLabel.Text = "Signed in as " + username;
                app.SignedInLabel.Visible = "on";

                app.SignOutButton.Visible = "on";

            else

                app.LoginFeedbackLabel.Text = "Wrong username or password.";

            end
        end

        % Button pushed function: SignOutButton
        function SignOutButtonPushed(app, event)
            appPath = phase_preprocessing_beta.projectRoot();

            resultFile = fullfile(appPath, "downloadasf", "login_result.json");
            requestFile = fullfile(appPath, "downloadasf", "login_request.json");

            if exist(resultFile, "file")
                delete(resultFile);
            end
            if exist(requestFile, "file")
                delete(requestFile);
            end

            app.EarthdataUsernameEditField.Value = "";
            app.EarthdataPasswordEditField.Value = "";

            app.EarthdataUsernameEditField.Visible = "on";
            app.EarthdataPasswordEditField.Visible = "on";

            app.LoginButton.Visible = "on";

            app.SignOutButton.Visible = "off";

            app.SignedInLabel.Visible = "off";



        end

        % Callback function
        function MapTypeDropDownValueChanged(app, event)

        end

        % Button pushed function: LoadLastDownloadButton
        function LoadLastDownloadButtonPushed(app, event)
            inputFile = fullfile(pwd, "downloadasf", "last_download_request.json");

            if ~exist(inputFile, "file")
                uialert(app.UIFigure, ...
                    "No previous download parameters found.", ...
                    "No Last Download", ...
                    "Icon", "warning");
                return;
            end

            params = jsondecode(fileread(inputFile));
            applyDownloaderParams(app, params);

            disp("Loaded last download parameters.");
        end

        % Button pushed function: TogglePreviewButton
        function TogglePreviewButtonPushed(app, event)
            if strcmp(app.PreviewPanel.Visible, "on")
                app.PreviewPanel.Visible = "off";
                app.TogglePreviewButton.Text = "Show Preview";
            else
                app.PreviewPanel.Visible = "on";
                app.TogglePreviewButton.Text = "Hide Preview";
            end
        end

        % Value changed function: SelectAllCheckBox
        function SelectAllCheckBoxValueChanged(app, event)
            tableData = app.PreviewResultsTable.Data;

            if isempty(tableData)
                return;
            end

            selectAll = app.SelectAllCheckBox.Value;

            for i = 1:size(tableData, 1)
                tableData{i,1} = selectAll;
            end

            app.PreviewResultsTable.Data = tableData;

            updatePreviewSelectionState(app);
            drawSelectedPreviewFootprints(app);
        end

        % Button pushed function: ToggleFootprintsButton
        function ToggleFootprintsButtonPushed(app, event)
            if isempty(app.SearchPreviewProducts)
                uialert(app.UIFigure, ...
                    "No footprints to show.", ...
                    "No Footprints", ...
                    "Icon", "info");
                return;
            end

            if app.FootprintsVisible

                clearPreviewFootprints(app);
                clearSelectedFootprints(app);

                app.FootprintsVisible = false;
                app.ToggleFootprintsButton.Text = "Show Footprints";

            else

                drawAllPreviewFootprints(app);
                drawSelectedPreviewFootprints(app);

                app.FootprintsVisible = true;
                app.ToggleFootprintsButton.Text = "Hide Footprints";

            end
        end

        % Cell edit callback: PreviewResultsTable
        function PreviewResultsTableCellSelection(app, event)
            % Only react when the Select checkbox column changes
            if event.Indices(2) ~= 1
                return;
            end

            updatePreviewSelectionState(app);
            drawSelectedPreviewFootprints(app);

        end

        % Button pushed function: DownloadSelectedButton
        function DownloadSelectedButtonPushed(app, event)
            if ~isEarthdataLoggedIn(app)
                uialert(app.UIFigure, ...
                    "Please log in before downloading from ASF.", ...
                    "Login Required", ...
                    "Icon", "warning");
                return;
            end

            disp("CALLING runBackendDownload")
            products = getSelectedPreviewProducts(app);

            if ~validateDownloadCompatibility(app, products)
                return;
            end

            if ~saveDownloadSelection(app, products)
                return;
            end

            if ~confirmDownloadSelection(app, products)
                return;
            end

            saveLastDownloadParams(app);

            totalBytes = sum([products.size_bytes]);

            runBackendDownload(app, totalBytes, products);
        end

        % Button pushed function: DownloadAllButton
        function DownloadAllButtonPushed(app, event)
            if ~isEarthdataLoggedIn(app)
                uialert(app.UIFigure, ...
                    "Please log in before downloading from ASF.", ...
                    "Login Required", ...
                    "Icon", "warning");
                return;
            end

            products = app.SearchPreviewProducts;

            if isempty(products)
                uialert(app.UIFigure, ...
                    "No images available to download.", ...
                    "Download Error", ...
                    "Icon", "warning");
                return;
            end

            if ~validateDownloadCompatibility(app, products)
                return;
            end

            if ~restoreAllDownloadFiles(app)
                return;
            end

            if ~confirmDownloadSelection(app, products)
                return;
            end

            saveLastDownloadParams(app);

            totalBytes = sum([products.size_bytes]);

            runBackendDownload(app, totalBytes, products);
        end

        % Button pushed function: ImportCSKButton
        function ImportCSKButtonPushed(app, event)
            % Select one or more Cosmo-SkyMed (.h5) images and copy them into the
            % project's "slaves" folder.
            [files, srcPath] = uigetfile('*.h5', ...
                'Select Cosmo-SkyMed images (.h5)', 'MultiSelect', 'on');
            if isequal(files, 0)
                return                      % user cancelled
            end
            if ischar(files)                % single selection -> cell
                files = {files};
            end

            % Destination: PHASE_Preprocessing\slaves
            currentFolder = phase_preprocessing_beta.projectRoot();
            slavesFolder = fullfile(currentFolder, 'PHASE_Preprocessing', 'slaves');
            if exist(slavesFolder, 'dir') ~= 7
                mkdir(slavesFolder);
            end

            % Copy with a progress dialog
            n = numel(files);
            nCopied = 0;
            dlg = uiprogressdlg(app.UIFigure, 'Title', 'Importing images', ...
                'Message', 'Copying...');
            for k = 1:n
                dlg.Value = k / n;
                dlg.Message = sprintf('Copying %d of %d: %s', k, n, files{k});
                src = fullfile(srcPath, files{k});
                dst = fullfile(slavesFolder, files{k});
                if exist(dst, 'file') == 2
                    continue                % already present -> skip
                end
                try
                    copyfile(src, dst);
                    nCopied = nCopied + 1;
                catch ME
                    uialert(app.UIFigure, ...
                        sprintf('Could not copy %s:\n%s', files{k}, ME.message), ...
                        'Copy error');
                end
            end
            close(dlg);

            refreshImportedImagesTable(app);

            app.initializeCSKMap(); % Refresh the AOI map with the new footprints

            updateOutput(app, sprintf('Imported %d new image(s) into slaves (%d selected).', nCopied, n));
        end

        % Button pushed function: OpenslavesfolderButton
        function OpenslavesfolderButtonPushed(app, event)
            currentFolder = phase_preprocessing_beta.projectRoot();
            slavesFolder = fullfile(currentFolder, 'PHASE_Preprocessing', 'slaves');
            if exist(slavesFolder, 'dir') ~= 7
                mkdir(slavesFolder);
            end
            winopen(slavesFolder);
        end

        % Value changed function: UpdateAlreadyProcessedDataCheckBox
        function UpdateAlreadyProcessedDataCheckBoxValueChanged(app, event)
            value = app.UpdateAlreadyProcessedDataCheckBox.Value;
            if value
                app.update_processed_data_SEN = 1;
                app.setSENUpdateControlsVisible(true);
                app.refreshSENUpdateContext();
            else
                app.update_processed_data_SEN = 0;
                app.setSENUpdateControlsVisible(false);
            end
        end

        % Value changed function: UpdateSearchEndDateDatePicker
        function UpdateSearchEndDateDatePickerValueChanged(app, event)
            if app.UpdateAlreadyProcessedDataCheckBox.Value
                app.update_search_results_SEN = struct([]);
                app.UpdateImagesTable.Data = {};
                app.DownloadUpdateImagesButton.Enable = 'off';
                app.DownloadRunUpdateImagesButton.Enable = 'off';
                app.UpdateSearchStatusLabel.Text = 'Click Search Images to query ASF.';
            end
        end

        % Button pushed function: SearchUpdateImagesButton
        function SearchUpdateImagesButtonPushed(app, event)
            app.searchSENUpdateImages();
        end

        % Button pushed function: DownloadUpdateImagesButton
        function DownloadUpdateImagesButtonPushed(app, event)
            app.downloadSelectedSENUpdateImages(false, true);
        end

        % Callback function
        function DownloadRunUpdateImagesDownloadRunUpdateImagesButtonPushed(app, event)
            app.downloadSelectedSENUpdateImages(true, false);
        end

        % Button pushed function: SelectAllUpdateImagesButton
        function SelectAllUpdateImagesButtonPushed(app, event)
            app.setAllSENUpdateImagesSelected(true);
        end

        % Button pushed function: DeselectAllUpdateImagesButton
        function DeselectAllUpdateImagesButtonPushed(app, event)
            app.setAllSENUpdateImagesSelected(false);
        end

        % Button pushed function: ImportSENButton
        function ImportSENButtonPushed(app, event)
            [files, pathName] = uigetfile({'*.zip', 'Sentinel-1 ZIP files (*.zip)'}, ...
                'Select Sentinel-1 images', 'MultiSelect', 'on');
            if isequal(files, 0)
                return;
            end
            if ischar(files)
                files = {files};
            end

            importMode = uiconfirm(app.UIFigure, ...
                ['Do you want to move or copy the selected Sentinel-1 ZIP files into PHASE_Preprocessing\slaves?' newline newline ...
                'Move is usually much faster when the files are on the same drive and does not duplicate disk usage.' newline ...
                'Copy keeps the original files, but requires additional disk space.'], ...
                'Import Sentinel-1 Images', ...
                'Options', {'Move', 'Copy', 'Cancel'}, ...
                'DefaultOption', 'Move', ...
                'CancelOption', 'Cancel', ...
                'Icon', 'question');
            if strcmp(importMode, 'Cancel')
                return;
            end

            currentFolder = phase_preprocessing_beta.projectRoot();
            slavesFolder = fullfile(currentFolder, 'PHASE_Preprocessing', 'slaves');
            if exist(slavesFolder, 'dir') ~= 7
                mkdir(slavesFolder);
            end

            useMove = strcmp(importMode, 'Move');
            if useMove
                progressMessage = 'Moving selected files into slaves...';
                pastVerb = 'Moved';
            else
                progressMessage = 'Copying selected files into slaves...';
                pastVerb = 'Copied';
            end

            importedCount = 0;
            skippedCount = 0;
            failedCount = 0;
            failedNames = {};
            dlg = uiprogressdlg(app.UIFigure, ...
                'Title', 'Importing Sentinel-1 images', ...
                'Message', progressMessage, ...
                'Indeterminate', 'off');
            for k = 1:numel(files)
                dlg.Value = k / numel(files);
                srcFile = fullfile(pathName, files{k});
                dstFile = fullfile(slavesFolder, files{k});
                if exist(dstFile, 'file') == 2
                    skippedCount = skippedCount + 1;
                    continue;
                end
                try
                    if useMove
                        [ok, msg] = movefile(srcFile, dstFile);
                    else
                        [ok, msg] = copyfile(srcFile, dstFile);
                    end
                    if ~ok
                        error('PHASE:ImportSENImageFailed', '%s', msg);
                    end
                    importedCount = importedCount + 1;
                catch ME
                    failedCount = failedCount + 1;
                    failedNames{end + 1} = sprintf('%s: %s', files{k}, ME.message); %#ok<AGROW>
                end
            end
            close(dlg);

            app.refreshImportedSENImagesTable();
            app.refreshSENUpdateContext();

            updateOutput(app, sprintf('%s %d Sentinel-1 image(s) into slaves. Skipped %d existing file(s). Failed %d file(s).', pastVerb, importedCount, skippedCount, failedCount));
            if failedCount > 0
                uialert(app.UIFigure, strjoin(failedNames, newline), 'Import Failed');
            end
        end

        % Button pushed function: OpenSENSlavesFolderButton
        function OpenSENSlavesFolderButtonPushed(app, event)
            currentFolder = phase_preprocessing_beta.projectRoot();
            slavesFolder = fullfile(currentFolder, 'PHASE_Preprocessing', 'slaves');
            if exist(slavesFolder, 'dir') ~= 7
                mkdir(slavesFolder);
            end

            if ispc
                winopen(slavesFolder);
            elseif ismac
                system(sprintf('open %s', app.quoteCommandArgument(slavesFolder)));
            else
                system(sprintf('xdg-open %s &', app.quoteCommandArgument(slavesFolder)));
            end
        end

        % Button pushed function: DownloadRunUpdateImagesButton
        function DownloadRunUpdateImagesButtonPushed(app, event)
            app.downloadSelectedSENUpdateImages(true, false);
        end
    end

    % Component initialization
    methods (Access = public)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = phase_preprocessing_beta.projectRoot();

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.Position = [100 100 1200 675];
            app.UIFigure.Name = 'PHASE · Preprocessing advanced workspace';

            % Create ConstellationSwitchLabel
            app.ConstellationSwitchLabel = uilabel(app.UIFigure);
            app.ConstellationSwitchLabel.HorizontalAlignment = 'center';
            app.ConstellationSwitchLabel.FontName = 'Manrope';
            app.ConstellationSwitchLabel.Position = [545 561 80 22];
            app.ConstellationSwitchLabel.Text = 'Constellation';

            % Create CosmoSkyMedPanel
            app.CosmoSkyMedPanel = uipanel(app.UIFigure);
            app.CosmoSkyMedPanel.Title = 'CosmoSkyMed';
            app.CosmoSkyMedPanel.Visible = 'off';
            app.CosmoSkyMedPanel.BackgroundColor = [1 1 1];
            app.CosmoSkyMedPanel.FontName = 'Manrope';
            app.CosmoSkyMedPanel.FontWeight = 'bold';
            app.CosmoSkyMedPanel.FontSize = 14;
            app.CosmoSkyMedPanel.Position = [1 1 1200 541];

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.CosmoSkyMedPanel);
            app.TabGroup2.SelectionChangedFcn = createCallbackFcn(app, @TabGroup2SelectionChanged, true);
            app.TabGroup2.Position = [0 0 1198 519];

            % Create GlobalVariablesTab_4
            app.GlobalVariablesTab_4 = uitab(app.TabGroup2);
            app.GlobalVariablesTab_4.Title = 'Global Variables';
            app.GlobalVariablesTab_4.BackgroundColor = [1 1 1];

            % Create CustomPythonEnvironmentLabel
            app.CustomPythonEnvironmentLabel = uilabel(app.GlobalVariablesTab_4);
            app.CustomPythonEnvironmentLabel.HorizontalAlignment = 'right';
            app.CustomPythonEnvironmentLabel.Position = [42 357 161 22];
            app.CustomPythonEnvironmentLabel.Text = 'Custom Python Environment:';

            % Create CustomPythonEnvironmentEditField_2
            app.CustomPythonEnvironmentEditField_2 = uieditfield(app.GlobalVariablesTab_4, 'text');
            app.CustomPythonEnvironmentEditField_2.ValueChangedFcn = createCallbackFcn(app, @CustomPythonEnvironmentEditField_2ValueChanged, true);
            app.CustomPythonEnvironmentEditField_2.FontName = 'Manrope';
            app.CustomPythonEnvironmentEditField_2.Position = [218 357 161 22];
            app.CustomPythonEnvironmentEditField_2.Value = 'python';

            % Create Label_4
            app.Label_4 = uilabel(app.GlobalVariablesTab_4);
            app.Label_4.HorizontalAlignment = 'center';
            app.Label_4.FontName = 'Manrope';
            app.Label_4.FontWeight = 'bold';
            app.Label_4.Position = [39 428 280 22];
            app.Label_4.Text = 'Name of the python 3.x environment in your os ';

            % Create PythonEnvironmentDropDown_2Label
            app.PythonEnvironmentDropDown_2Label = uilabel(app.GlobalVariablesTab_4);
            app.PythonEnvironmentDropDown_2Label.HorizontalAlignment = 'right';
            app.PythonEnvironmentDropDown_2Label.Position = [44 397 113 22];
            app.PythonEnvironmentDropDown_2Label.Text = 'Python Environment';

            % Create PythonEnvironmentDropDown_2
            app.PythonEnvironmentDropDown_2 = uidropdown(app.GlobalVariablesTab_4);
            app.PythonEnvironmentDropDown_2.Items = {'python', 'python3', 'python3.11', 'Other'};
            app.PythonEnvironmentDropDown_2.ValueChangedFcn = createCallbackFcn(app, @PythonEnvironmentDropDown_2ValueChanged, true);
            app.PythonEnvironmentDropDown_2.FontName = 'Manrope';
            app.PythonEnvironmentDropDown_2.BackgroundColor = [1 1 1];
            app.PythonEnvironmentDropDown_2.Position = [172 397 100 22];
            app.PythonEnvironmentDropDown_2.Value = 'python';

            % Create AOITab
            app.AOITab = uitab(app.TabGroup2);
            app.AOITab.Title = 'AOI';
            app.AOITab.BackgroundColor = [1 1 1];

            % Create Panel
            app.Panel = uipanel(app.AOITab);
            app.Panel.BorderType = 'none';
            app.Panel.BackgroundColor = [1 1 1];
            app.Panel.FontName = 'Manrope';
            app.Panel.Scrollable = 'on';
            app.Panel.Position = [34 32 834 436];

            % Create DrawCSKAOIButton
            app.DrawCSKAOIButton = uibutton(app.AOITab, 'push');
            app.DrawCSKAOIButton.ButtonPushedFcn = createCallbackFcn(app, @DrawCSKAOIButtonPushed, true);
            app.DrawCSKAOIButton.BackgroundColor = [1 1 1];
            app.DrawCSKAOIButton.FontName = 'Manrope';
            app.DrawCSKAOIButton.Position = [899 445 149 23];
            app.DrawCSKAOIButton.Text = 'Reset Map';

            % Create MinlongitudeEditField_2Label
            app.MinlongitudeEditField_2Label = uilabel(app.AOITab);
            app.MinlongitudeEditField_2Label.HorizontalAlignment = 'right';
            app.MinlongitudeEditField_2Label.Position = [899 131 76 22];
            app.MinlongitudeEditField_2Label.Text = 'Min longitude';

            % Create MinlongitudeEditField_2
            app.MinlongitudeEditField_2 = uieditfield(app.AOITab, 'numeric');
            app.MinlongitudeEditField_2.Limits = [-180 180];
            app.MinlongitudeEditField_2.ValueDisplayFormat = '%.3f';
            app.MinlongitudeEditField_2.ValueChangedFcn = createCallbackFcn(app, @MinlongitudeEditField_2ValueChanged, true);
            app.MinlongitudeEditField_2.Editable = 'off';
            app.MinlongitudeEditField_2.FontName = 'Manrope';
            app.MinlongitudeEditField_2.Position = [994 131 100 22];
            app.MinlongitudeEditField_2.Value = -180;

            % Create MaxlongitudeEditField_2Label
            app.MaxlongitudeEditField_2Label = uilabel(app.AOITab);
            app.MaxlongitudeEditField_2Label.HorizontalAlignment = 'right';
            app.MaxlongitudeEditField_2Label.Position = [899 98 80 22];
            app.MaxlongitudeEditField_2Label.Text = 'Max longitude';

            % Create MaxlongitudeEditField_2
            app.MaxlongitudeEditField_2 = uieditfield(app.AOITab, 'numeric');
            app.MaxlongitudeEditField_2.Limits = [-180 180];
            app.MaxlongitudeEditField_2.ValueDisplayFormat = '%.3f';
            app.MaxlongitudeEditField_2.ValueChangedFcn = createCallbackFcn(app, @MaxlongitudeEditField_2ValueChanged, true);
            app.MaxlongitudeEditField_2.Editable = 'off';
            app.MaxlongitudeEditField_2.FontName = 'Manrope';
            app.MaxlongitudeEditField_2.Position = [994 98 100 22];
            app.MaxlongitudeEditField_2.Value = 180;

            % Create MaxlatitudeEditField_2Label
            app.MaxlatitudeEditField_2Label = uilabel(app.AOITab);
            app.MaxlatitudeEditField_2Label.HorizontalAlignment = 'right';
            app.MaxlatitudeEditField_2Label.Position = [899 32 70 22];
            app.MaxlatitudeEditField_2Label.Text = 'Max latitude';

            % Create MaxlatitudeEditField_2
            app.MaxlatitudeEditField_2 = uieditfield(app.AOITab, 'numeric');
            app.MaxlatitudeEditField_2.Limits = [-90 90];
            app.MaxlatitudeEditField_2.ValueDisplayFormat = '%.3f';
            app.MaxlatitudeEditField_2.ValueChangedFcn = createCallbackFcn(app, @MaxlatitudeEditField_2ValueChanged, true);
            app.MaxlatitudeEditField_2.Editable = 'off';
            app.MaxlatitudeEditField_2.FontName = 'Manrope';
            app.MaxlatitudeEditField_2.Position = [994 32 100 22];
            app.MaxlatitudeEditField_2.Value = 90;

            % Create MinlatitudeEditField_2Label
            app.MinlatitudeEditField_2Label = uilabel(app.AOITab);
            app.MinlatitudeEditField_2Label.HorizontalAlignment = 'right';
            app.MinlatitudeEditField_2Label.Position = [899 64 66 22];
            app.MinlatitudeEditField_2Label.Text = 'Min latitude';

            % Create MinlatitudeEditField_2
            app.MinlatitudeEditField_2 = uieditfield(app.AOITab, 'numeric');
            app.MinlatitudeEditField_2.Limits = [-90 90];
            app.MinlatitudeEditField_2.ValueDisplayFormat = '%.3f';
            app.MinlatitudeEditField_2.ValueChangedFcn = createCallbackFcn(app, @MinlatitudeEditField_2ValueChanged, true);
            app.MinlatitudeEditField_2.Editable = 'off';
            app.MinlatitudeEditField_2.FontName = 'Manrope';
            app.MinlatitudeEditField_2.Position = [994 64 100 22];
            app.MinlatitudeEditField_2.Value = -90;

            % Create AOIboxboundariesLabel_2
            app.AOIboxboundariesLabel_2 = uilabel(app.AOITab);
            app.AOIboxboundariesLabel_2.HorizontalAlignment = 'center';
            app.AOIboxboundariesLabel_2.FontName = 'Manrope';
            app.AOIboxboundariesLabel_2.FontWeight = 'bold';
            app.AOIboxboundariesLabel_2.Position = [899 167 119 22];
            app.AOIboxboundariesLabel_2.Text = 'AOI box boundaries';

            % Create MasterProcessingTab_2
            app.MasterProcessingTab_2 = uitab(app.TabGroup2);
            app.MasterProcessingTab_2.Title = 'Master Processing';
            app.MasterProcessingTab_2.BackgroundColor = [1 1 1];

            % Create MasterdateDatePicker_2Label
            app.MasterdateDatePicker_2Label = uilabel(app.MasterProcessingTab_2);
            app.MasterdateDatePicker_2Label.HorizontalAlignment = 'right';
            app.MasterdateDatePicker_2Label.FontWeight = 'bold';
            app.MasterdateDatePicker_2Label.Position = [40 437 72 22];
            app.MasterdateDatePicker_2Label.Text = 'Master date';

            % Create MasterdateDatePicker_2
            app.MasterdateDatePicker_2 = uidatepicker(app.MasterProcessingTab_2);
            app.MasterdateDatePicker_2.ValueChangedFcn = createCallbackFcn(app, @MasterdateDatePicker_2ValueChanged, true);
            app.MasterdateDatePicker_2.FontName = 'Manrope';
            app.MasterdateDatePicker_2.Enable = 'off';
            app.MasterdateDatePicker_2.Position = [127 437 150 22];
            app.MasterdateDatePicker_2.Value = datetime([2020 3 5]);

            % Create MasterprocessingCheckBox_2
            app.MasterprocessingCheckBox_2 = uicheckbox(app.MasterProcessingTab_2);
            app.MasterprocessingCheckBox_2.ValueChangedFcn = createCallbackFcn(app, @MasterprocessingCheckBox_2ValueChanged, true);
            app.MasterprocessingCheckBox_2.Text = 'Master processing';
            app.MasterprocessingCheckBox_2.FontName = 'Manrope';
            app.MasterprocessingCheckBox_2.FontWeight = 'bold';
            app.MasterprocessingCheckBox_2.Position = [45 386 130 22];
            app.MasterprocessingCheckBox_2.Value = true;

            % Create AutoMasterCheckBox_2
            app.AutoMasterCheckBox_2 = uicheckbox(app.MasterProcessingTab_2);
            app.AutoMasterCheckBox_2.ValueChangedFcn = createCallbackFcn(app, @AutoMasterCheckBox_2ValueChanged, true);
            app.AutoMasterCheckBox_2.Text = 'Auto-select master';
            app.AutoMasterCheckBox_2.FontName = 'Manrope';
            app.AutoMasterCheckBox_2.Position = [337 435 129 22];
            app.AutoMasterCheckBox_2.Value = true;

            % Create SlavesProcessingTab_2
            app.SlavesProcessingTab_2 = uitab(app.TabGroup2);
            app.SlavesProcessingTab_2.Title = 'Slaves Processing';
            app.SlavesProcessingTab_2.BackgroundColor = [1 1 1];

            % Create SlavesremovalafterprocessingCheckBox_2
            app.SlavesremovalafterprocessingCheckBox_2 = uicheckbox(app.SlavesProcessingTab_2);
            app.SlavesremovalafterprocessingCheckBox_2.ValueChangedFcn = createCallbackFcn(app, @SlavesremovalafterprocessingCheckBox_2ValueChanged, true);
            app.SlavesremovalafterprocessingCheckBox_2.Text = 'Slaves removal after processing';
            app.SlavesremovalafterprocessingCheckBox_2.FontName = 'Manrope';
            app.SlavesremovalafterprocessingCheckBox_2.FontWeight = 'bold';
            app.SlavesremovalafterprocessingCheckBox_2.Position = [48 437 207 22];

            % Create DEMinterferogramDropDown_2Label
            app.DEMinterferogramDropDown_2Label = uilabel(app.SlavesProcessingTab_2);
            app.DEMinterferogramDropDown_2Label.HorizontalAlignment = 'right';
            app.DEMinterferogramDropDown_2Label.FontWeight = 'bold';
            app.DEMinterferogramDropDown_2Label.Position = [41 309 113 22];
            app.DEMinterferogramDropDown_2Label.Text = 'DEM interferogram';

            % Create DEMinterferogramDropDown_2
            app.DEMinterferogramDropDown_2 = uidropdown(app.SlavesProcessingTab_2);
            app.DEMinterferogramDropDown_2.Items = {'SRTM 1Sec HGT', 'SRTM 3Sec', 'Copernicus 30m Global DEM', 'Copernicus 90m Global DEM', 'CDEM', 'GETASSE30', 'External DEM'};
            app.DEMinterferogramDropDown_2.ValueChangedFcn = createCallbackFcn(app, @DEMinterferogramDropDown_2ValueChanged, true);
            app.DEMinterferogramDropDown_2.FontName = 'Manrope';
            app.DEMinterferogramDropDown_2.FontWeight = 'bold';
            app.DEMinterferogramDropDown_2.BackgroundColor = [1 1 1];
            app.DEMinterferogramDropDown_2.Position = [169 309 198 22];
            app.DEMinterferogramDropDown_2.Value = 'SRTM 1Sec HGT';

            % Create DEMifgpathEditField_2Label
            app.DEMifgpathEditField_2Label = uilabel(app.SlavesProcessingTab_2);
            app.DEMifgpathEditField_2Label.HorizontalAlignment = 'right';
            app.DEMifgpathEditField_2Label.Position = [52 267 74 22];
            app.DEMifgpathEditField_2Label.Text = 'DEM ifg path';

            % Create DEMifgpathEditField_2
            app.DEMifgpathEditField_2 = uieditfield(app.SlavesProcessingTab_2, 'text');
            app.DEMifgpathEditField_2.ValueChangedFcn = createCallbackFcn(app, @DEMifgpathEditField_2ValueChanged, true);
            app.DEMifgpathEditField_2.HorizontalAlignment = 'right';
            app.DEMifgpathEditField_2.FontName = 'Manrope';
            app.DEMifgpathEditField_2.Position = [141 267 502 22];

            % Create onlyforExternalDEMLabel_3
            app.onlyforExternalDEMLabel_3 = uilabel(app.SlavesProcessingTab_2);
            app.onlyforExternalDEMLabel_3.FontName = 'Manrope';
            app.onlyforExternalDEMLabel_3.Position = [779 263 131 22];
            app.onlyforExternalDEMLabel_3.Text = '(only for External DEM)';

            % Create FirststepDropDown_2Label
            app.FirststepDropDown_2Label = uilabel(app.SlavesProcessingTab_2);
            app.FirststepDropDown_2Label.HorizontalAlignment = 'right';
            app.FirststepDropDown_2Label.FontWeight = 'bold';
            app.FirststepDropDown_2Label.Position = [41 396 59 22];
            app.FirststepDropDown_2Label.Text = 'First step';

            % Create FirststepDropDown_2
            app.FirststepDropDown_2 = uidropdown(app.SlavesProcessingTab_2);
            app.FirststepDropDown_2.Items = {'1', '2', '3', '4', '5', '6'};
            app.FirststepDropDown_2.ValueChangedFcn = createCallbackFcn(app, @FirststepDropDown_2ValueChanged, true);
            app.FirststepDropDown_2.FontName = 'Manrope';
            app.FirststepDropDown_2.FontWeight = 'bold';
            app.FirststepDropDown_2.BackgroundColor = [1 1 1];
            app.FirststepDropDown_2.Position = [115 396 62 22];
            app.FirststepDropDown_2.Value = '1';

            % Create CoregistrationGCPsnumberEditFieldLabel
            app.CoregistrationGCPsnumberEditFieldLabel = uilabel(app.SlavesProcessingTab_2);
            app.CoregistrationGCPsnumberEditFieldLabel.HorizontalAlignment = 'right';
            app.CoregistrationGCPsnumberEditFieldLabel.FontWeight = 'bold';
            app.CoregistrationGCPsnumberEditFieldLabel.Position = [41 355 170 22];
            app.CoregistrationGCPsnumberEditFieldLabel.Text = 'Coregistration GCPs number';

            % Create CoregistrationGCPsnumberEditField
            app.CoregistrationGCPsnumberEditField = uieditfield(app.SlavesProcessingTab_2, 'numeric');
            app.CoregistrationGCPsnumberEditField.Limits = [0 Inf];
            app.CoregistrationGCPsnumberEditField.ValueDisplayFormat = '%.0f';
            app.CoregistrationGCPsnumberEditField.ValueChangedFcn = createCallbackFcn(app, @CoregistrationGCPsnumberEditFieldValueChanged, true);
            app.CoregistrationGCPsnumberEditField.FontName = 'Manrope';
            app.CoregistrationGCPsnumberEditField.FontWeight = 'bold';
            app.CoregistrationGCPsnumberEditField.Position = [226 355 100 22];
            app.CoregistrationGCPsnumberEditField.Value = 10000;

            % Create BrowseButton_3
            app.BrowseButton_3 = uibutton(app.SlavesProcessingTab_2, 'push');
            app.BrowseButton_3.ButtonPushedFcn = createCallbackFcn(app, @BrowseButton_3Pushed, true);
            app.BrowseButton_3.BackgroundColor = [1 1 1];
            app.BrowseButton_3.FontName = 'Manrope';
            app.BrowseButton_3.Position = [662 266 100 23];
            app.BrowseButton_3.Text = 'Browse';

            % Create CoherenceandLIATab_2
            app.CoherenceandLIATab_2 = uitab(app.TabGroup2);
            app.CoherenceandLIATab_2.Title = 'Coherence and LIA';
            app.CoherenceandLIATab_2.BackgroundColor = [1 1 1];

            % Create TerraincorrectedCoherenceandLIACheckBox_2
            app.TerraincorrectedCoherenceandLIACheckBox_2 = uicheckbox(app.CoherenceandLIATab_2);
            app.TerraincorrectedCoherenceandLIACheckBox_2.ValueChangedFcn = createCallbackFcn(app, @TerraincorrectedCoherenceandLIACheckBox_2ValueChanged, true);
            app.TerraincorrectedCoherenceandLIACheckBox_2.Text = 'Terrain-corrected Coherence and LIA';
            app.TerraincorrectedCoherenceandLIACheckBox_2.FontName = 'Manrope';
            app.TerraincorrectedCoherenceandLIACheckBox_2.FontWeight = 'bold';
            app.TerraincorrectedCoherenceandLIACheckBox_2.Position = [48 437 238 22];
            app.TerraincorrectedCoherenceandLIACheckBox_2.Value = true;

            % Create EPSGcodeEditField_2Label
            app.EPSGcodeEditField_2Label = uilabel(app.CoherenceandLIATab_2);
            app.EPSGcodeEditField_2Label.HorizontalAlignment = 'right';
            app.EPSGcodeEditField_2Label.Position = [52 364 68 22];
            app.EPSGcodeEditField_2Label.Text = 'EPSG code';

            % Create EPSGcodeEditField_2
            app.EPSGcodeEditField_2 = uieditfield(app.CoherenceandLIATab_2, 'numeric');
            app.EPSGcodeEditField_2.Limits = [0 99999];
            app.EPSGcodeEditField_2.ValueDisplayFormat = '%.0f';
            app.EPSGcodeEditField_2.ValueChangedFcn = createCallbackFcn(app, @EPSGcodeEditField_2ValueChanged, true);
            app.EPSGcodeEditField_2.FontName = 'Manrope';
            app.EPSGcodeEditField_2.Position = [135 364 100 22];
            app.EPSGcodeEditField_2.Value = 32633;

            % Create Label_5
            app.Label_5 = uilabel(app.CoherenceandLIATab_2);
            app.Label_5.FontName = 'Manrope';
            app.Label_5.FontWeight = 'bold';
            app.Label_5.Position = [48 395 670 22];
            app.Label_5.Text = 'EPSG code of the cartographic projection of the terrain-corrected .TIFF files (choose accordingly to the location)';

            % Create ComputationalResourcesTab_2
            app.ComputationalResourcesTab_2 = uitab(app.TabGroup2);
            app.ComputationalResourcesTab_2.Title = 'Computational Resources';
            app.ComputationalResourcesTab_2.BackgroundColor = [1 1 1];

            % Create FullpathtoSNAPgptfolderLabel_2
            app.FullpathtoSNAPgptfolderLabel_2 = uilabel(app.ComputationalResourcesTab_2);
            app.FullpathtoSNAPgptfolderLabel_2.FontName = 'Manrope';
            app.FullpathtoSNAPgptfolderLabel_2.FontWeight = 'bold';
            app.FullpathtoSNAPgptfolderLabel_2.Position = [48 437 165 22];
            app.FullpathtoSNAPgptfolderLabel_2.Text = 'Full path to SNAP gpt folder';

            % Create PathEditField_2Label
            app.PathEditField_2Label = uilabel(app.ComputationalResourcesTab_2);
            app.PathEditField_2Label.HorizontalAlignment = 'right';
            app.PathEditField_2Label.Position = [54 406 30 22];
            app.PathEditField_2Label.Text = 'Path';

            % Create PathEditField_2
            app.PathEditField_2 = uieditfield(app.ComputationalResourcesTab_2, 'text');
            app.PathEditField_2.ValueChangedFcn = createCallbackFcn(app, @PathEditField_2ValueChanged, true);
            app.PathEditField_2.FontName = 'Manrope';
            app.PathEditField_2.Position = [99 406 284 22];
            app.PathEditField_2.Value = 'C:\Program Files\snap\bin\gpt';

            % Create NumberofcorestobeusedintheprocessingLabel_2
            app.NumberofcorestobeusedintheprocessingLabel_2 = uilabel(app.ComputationalResourcesTab_2);
            app.NumberofcorestobeusedintheprocessingLabel_2.FontName = 'Manrope';
            app.NumberofcorestobeusedintheprocessingLabel_2.FontWeight = 'bold';
            app.NumberofcorestobeusedintheprocessingLabel_2.Position = [45 359 269 22];
            app.NumberofcorestobeusedintheprocessingLabel_2.Text = 'Number of cores to be used in the processing';

            % Create CPUEditField_2Label
            app.CPUEditField_2Label = uilabel(app.ComputationalResourcesTab_2);
            app.CPUEditField_2Label.HorizontalAlignment = 'right';
            app.CPUEditField_2Label.Position = [57 325 30 22];
            app.CPUEditField_2Label.Text = 'CPU';

            % Create CPUEditField_2
            app.CPUEditField_2 = uieditfield(app.ComputationalResourcesTab_2, 'numeric');
            app.CPUEditField_2.Limits = [0 1000];
            app.CPUEditField_2.RoundFractionalValues = 'on';
            app.CPUEditField_2.ValueDisplayFormat = '%.0f';
            app.CPUEditField_2.ValueChangedFcn = createCallbackFcn(app, @CPUEditField_2ValueChanged, true);
            app.CPUEditField_2.HorizontalAlignment = 'left';
            app.CPUEditField_2.FontName = 'Manrope';
            app.CPUEditField_2.Position = [102 325 73 22];
            app.CPUEditField_2.Value = 8;

            % Create RAMtobeusedintheprocessingGBformatnnGLabel_2
            app.RAMtobeusedintheprocessingGBformatnnGLabel_2 = uilabel(app.ComputationalResourcesTab_2);
            app.RAMtobeusedintheprocessingGBformatnnGLabel_2.FontName = 'Manrope';
            app.RAMtobeusedintheprocessingGBformatnnGLabel_2.FontWeight = 'bold';
            app.RAMtobeusedintheprocessingGBformatnnGLabel_2.Position = [48 273 320 22];
            app.RAMtobeusedintheprocessingGBformatnnGLabel_2.Text = 'RAM to be used in the processing [GB] (format "nnG")';

            % Create CacheEditField_2Label
            app.CacheEditField_2Label = uilabel(app.ComputationalResourcesTab_2);
            app.CacheEditField_2Label.HorizontalAlignment = 'right';
            app.CacheEditField_2Label.Position = [55 242 40 22];
            app.CacheEditField_2Label.Text = 'Cache';

            % Create CacheEditField_2
            app.CacheEditField_2 = uieditfield(app.ComputationalResourcesTab_2, 'text');
            app.CacheEditField_2.ValueChangedFcn = createCallbackFcn(app, @CacheEditField_2ValueChanged, true);
            app.CacheEditField_2.FontName = 'Manrope';
            app.CacheEditField_2.Position = [110 242 100 22];
            app.CacheEditField_2.Value = '26G';

            % Create ImagesTab_CSK
            app.ImagesTab_CSK = uitab(app.TabGroup2);
            app.ImagesTab_CSK.Title = 'Images';
            app.ImagesTab_CSK.BackgroundColor = [1 1 1];

            % Create ImportedImagesTable
            app.ImportedImagesTable = uitable(app.ImagesTab_CSK);
            app.ImportedImagesTable.ColumnName = {'File'; 'Date'; 'Type'; 'Status'};
            app.ImportedImagesTable.RowName = {};
            app.ImportedImagesTable.Position = [261 64 675 325];

            % Create Label_10
            app.Label_10 = uilabel(app.ImagesTab_CSK);
            app.Label_10.HorizontalAlignment = 'center';
            app.Label_10.FontName = 'Manrope';
            app.Label_10.FontSize = 10;
            app.Label_10.Position = [425 402 347 22];
            app.Label_10.Text = 'Selected images will be copied into the project''s slaves folder';

            % Create ImportCSKButton
            app.ImportCSKButton = uibutton(app.ImagesTab_CSK, 'push');
            app.ImportCSKButton.ButtonPushedFcn = createCallbackFcn(app, @ImportCSKButtonPushed, true);
            app.ImportCSKButton.FontName = 'Manrope';
            app.ImportCSKButton.FontSize = 14;
            app.ImportCSKButton.Position = [399 431 399 39];
            app.ImportCSKButton.Text = 'Import COSMO-SkyMed images (.h5)';

            % Create OpenslavesfolderButton
            app.OpenslavesfolderButton = uibutton(app.ImagesTab_CSK, 'push');
            app.OpenslavesfolderButton.ButtonPushedFcn = createCallbackFcn(app, @OpenslavesfolderButtonPushed, true);
            app.OpenslavesfolderButton.FontName = 'Manrope';
            app.OpenslavesfolderButton.Position = [521 20 155 24];
            app.OpenslavesfolderButton.Text = 'Open slaves folder';

            % Create SaveLoadTab_2
            app.SaveLoadTab_2 = uitab(app.TabGroup2);
            app.SaveLoadTab_2.Title = 'Save/Load';
            app.SaveLoadTab_2.BackgroundColor = [1 1 1];

            % Create SaveButton_2
            app.SaveButton_2 = uibutton(app.SaveLoadTab_2, 'push');
            app.SaveButton_2.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton_2.BackgroundColor = [1 1 1];
            app.SaveButton_2.FontName = 'Manrope';
            app.SaveButton_2.Position = [44 404 100 23];
            app.SaveButton_2.Text = 'Save';

            % Create SavetheconfiguredparametersfortheInSARpreprocessingLabel_2
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_2 = uilabel(app.SaveLoadTab_2);
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_2.FontName = 'Manrope';
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_2.FontWeight = 'bold';
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_2.Position = [44 439 367 22];
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_2.Text = 'Save the configured parameters for the InSAR pre-processing';

            % Create StatusLamp_2Label
            app.StatusLamp_2Label = uilabel(app.SaveLoadTab_2);
            app.StatusLamp_2Label.HorizontalAlignment = 'right';
            app.StatusLamp_2Label.Position = [171 403 39 22];
            app.StatusLamp_2Label.Text = 'Status';

            % Create StatusLamp_2
            app.StatusLamp_2 = uilamp(app.SaveLoadTab_2);
            app.StatusLamp_2.Position = [225 403 20 20];
            app.StatusLamp_2.Color = [1 0 0];

            % Create LoadButton_2
            app.LoadButton_2 = uibutton(app.SaveLoadTab_2, 'push');
            app.LoadButton_2.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton_2.BackgroundColor = [1 1 1];
            app.LoadButton_2.FontName = 'Manrope';
            app.LoadButton_2.Position = [42 291 100 23];
            app.LoadButton_2.Text = 'Load';

            % Create SavetheconfiguredparametersfortheInSARpreprocessingLabel_4
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_4 = uilabel(app.SaveLoadTab_2);
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_4.FontName = 'Manrope';
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_4.FontWeight = 'bold';
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_4.Position = [44 326 502 22];
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_4.Text = 'Load the parameters for the InSAR pre-processing from the already existing .mat file';

            % Create StatusLamp_4Label
            app.StatusLamp_4Label = uilabel(app.SaveLoadTab_2);
            app.StatusLamp_4Label.HorizontalAlignment = 'right';
            app.StatusLamp_4Label.Position = [170 292 39 22];
            app.StatusLamp_4Label.Text = 'Status';

            % Create StatusLamp_4
            app.StatusLamp_4 = uilamp(app.SaveLoadTab_2);
            app.StatusLamp_4.Position = [224 292 20 20];
            app.StatusLamp_4.Color = [1 0 0];

            % Create RunTab_2
            app.RunTab_2 = uitab(app.TabGroup2);
            app.RunTab_2.Title = 'Run';
            app.RunTab_2.BackgroundColor = [1 1 1];

            % Create Label_8
            app.Label_8 = uilabel(app.RunTab_2);
            app.Label_8.FontName = 'Manrope';
            app.Label_8.FontWeight = 'bold';
            app.Label_8.Position = [40 439 525 22];
            app.Label_8.Text = 'Press the button below to start the pre-processing of the SAR images for the PS analysis';

            % Create StartButton_2
            app.StartButton_2 = uibutton(app.RunTab_2, 'push');
            app.StartButton_2.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton_2.BackgroundColor = [0.9412 0.9412 0.9412];
            app.StartButton_2.FontName = 'Manrope';
            app.StartButton_2.Position = [40 406 100 23];
            app.StartButton_2.Text = 'Start';

            % Create Label_9
            app.Label_9 = uilabel(app.RunTab_2);
            app.Label_9.FontName = 'Manrope';
            app.Label_9.FontWeight = 'bold';
            app.Label_9.Position = [40 342 117 22];
            app.Label_9.Text = 'Execution outputs:';

            % Create PreprocessingstatusLamp_2Label
            app.PreprocessingstatusLamp_2Label = uilabel(app.RunTab_2);
            app.PreprocessingstatusLamp_2Label.HorizontalAlignment = 'right';
            app.PreprocessingstatusLamp_2Label.Position = [181 405 117 22];
            app.PreprocessingstatusLamp_2Label.Text = 'Preprocessing status';

            % Create PreprocessingstatusLamp_2
            app.PreprocessingstatusLamp_2 = uilamp(app.RunTab_2);
            app.PreprocessingstatusLamp_2.Position = [313 405 20 20];
            app.PreprocessingstatusLamp_2.Color = [1 0 0];

            % Create MessagesTextAreaLabel
            app.MessagesTextAreaLabel = uilabel(app.RunTab_2);
            app.MessagesTextAreaLabel.HorizontalAlignment = 'right';
            app.MessagesTextAreaLabel.Position = [34 307 60 22];
            app.MessagesTextAreaLabel.Text = 'Messages';

            % Create MessagesTextArea_2
            app.MessagesTextArea_2 = uitextarea(app.RunTab_2);
            app.MessagesTextArea_2.Editable = 'off';
            app.MessagesTextArea_2.FontName = 'Manrope';
            app.MessagesTextArea_2.Position = [109 91 495 240];

            % Create StopButton_2
            app.StopButton_2 = uibutton(app.RunTab_2, 'push');
            app.StopButton_2.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);
            app.StopButton_2.BackgroundColor = [1 0.8588 0.8588];
            app.StopButton_2.FontName = 'Manrope';
            app.StopButton_2.Position = [631 309 100 23];
            app.StopButton_2.Text = 'Stop';

            % Create ThecodewillstopattheendofthecurrentstepLabel
            app.ThecodewillstopattheendofthecurrentstepLabel = uilabel(app.RunTab_2);
            app.ThecodewillstopattheendofthecurrentstepLabel.FontName = 'Manrope';
            app.ThecodewillstopattheendofthecurrentstepLabel.Position = [749 309 271 22];
            app.ThecodewillstopattheendofthecurrentstepLabel.Text = 'The code will stop at the end of the current step';

            % Create Sentinel1Panel
            app.Sentinel1Panel = uipanel(app.UIFigure);
            app.Sentinel1Panel.Title = 'Sentinel1';
            app.Sentinel1Panel.Visible = 'off';
            app.Sentinel1Panel.BackgroundColor = [1 1 1];
            app.Sentinel1Panel.FontName = 'Manrope';
            app.Sentinel1Panel.FontWeight = 'bold';
            app.Sentinel1Panel.FontSize = 14;
            app.Sentinel1Panel.Position = [1 1 1200 541];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.Sentinel1Panel);
            app.TabGroup.SelectionChangedFcn = createCallbackFcn(app, @TabGroupSelectionChanged, true);
            app.TabGroup.Position = [0 -1 1199 520];

            % Create DownloadTab
            app.DownloadTab = uitab(app.TabGroup);
            app.DownloadTab.Title = 'Download';
            app.DownloadTab.BackgroundColor = [1 1 1];

            % Create DownloaderMapPanel
            app.DownloaderMapPanel = uipanel(app.DownloadTab);
            app.DownloaderMapPanel.BorderType = 'none';
            app.DownloaderMapPanel.BackgroundColor = [1 1 1];
            app.DownloaderMapPanel.Position = [-67 -4 960 464];

            % Create FilterPanel
            app.FilterPanel = uipanel(app.DownloadTab);
            app.FilterPanel.Title = 'Filter';
            app.FilterPanel.BackgroundColor = [1 1 1];
            app.FilterPanel.Position = [5 54 792 390];

            % Create PathStartEditFieldLabel
            app.PathStartEditFieldLabel = uilabel(app.FilterPanel);
            app.PathStartEditFieldLabel.HorizontalAlignment = 'right';
            app.PathStartEditFieldLabel.FontColor = [0.149 0.149 0.149];
            app.PathStartEditFieldLabel.Position = [367 156 58 22];
            app.PathStartEditFieldLabel.Text = 'Path Start';

            % Create PathStartEditField
            app.PathStartEditField = uieditfield(app.FilterPanel, 'numeric');
            app.PathStartEditField.AllowEmpty = 'on';
            app.PathStartEditField.FontColor = [0.149 0.149 0.149];
            app.PathStartEditField.Position = [440 156 94 22];
            app.PathStartEditField.Value = [];

            % Create PathEndEditFieldLabel
            app.PathEndEditFieldLabel = uilabel(app.FilterPanel);
            app.PathEndEditFieldLabel.HorizontalAlignment = 'right';
            app.PathEndEditFieldLabel.FontColor = [0.149 0.149 0.149];
            app.PathEndEditFieldLabel.Position = [557 156 54 22];
            app.PathEndEditFieldLabel.Text = 'Path End';

            % Create PathEndEditField
            app.PathEndEditField = uieditfield(app.FilterPanel, 'numeric');
            app.PathEndEditField.AllowEmpty = 'on';
            app.PathEndEditField.FontColor = [0.149 0.149 0.149];
            app.PathEndEditField.Position = [626 156 96 22];
            app.PathEndEditField.Value = [];

            % Create FrameStartEditFieldLabel
            app.FrameStartEditFieldLabel = uilabel(app.FilterPanel);
            app.FrameStartEditFieldLabel.HorizontalAlignment = 'right';
            app.FrameStartEditFieldLabel.FontColor = [0.149 0.149 0.149];
            app.FrameStartEditFieldLabel.Position = [367 126 68 22];
            app.FrameStartEditFieldLabel.Text = 'Frame Start';

            % Create FrameStartEditField
            app.FrameStartEditField = uieditfield(app.FilterPanel, 'numeric');
            app.FrameStartEditField.AllowEmpty = 'on';
            app.FrameStartEditField.FontColor = [0.149 0.149 0.149];
            app.FrameStartEditField.Position = [450 126 84 22];
            app.FrameStartEditField.Value = [];

            % Create FrameEndEditFieldLabel
            app.FrameEndEditFieldLabel = uilabel(app.FilterPanel);
            app.FrameEndEditFieldLabel.HorizontalAlignment = 'right';
            app.FrameEndEditFieldLabel.FontColor = [0.149 0.149 0.149];
            app.FrameEndEditFieldLabel.Position = [557 126 64 22];
            app.FrameEndEditFieldLabel.Text = 'Frame End';

            % Create FrameEndEditField
            app.FrameEndEditField = uieditfield(app.FilterPanel, 'numeric');
            app.FrameEndEditField.AllowEmpty = 'on';
            app.FrameEndEditField.FontColor = [0.149 0.149 0.149];
            app.FrameEndEditField.Position = [636 126 86 22];
            app.FrameEndEditField.Value = [];

            % Create DatasetDropDownLabel
            app.DatasetDropDownLabel = uilabel(app.FilterPanel);
            app.DatasetDropDownLabel.HorizontalAlignment = 'right';
            app.DatasetDropDownLabel.Position = [10 336 46 22];
            app.DatasetDropDownLabel.Text = 'Dataset';

            % Create DatasetDropDown
            app.DatasetDropDown = uidropdown(app.FilterPanel);
            app.DatasetDropDown.Items = {'SENTINEL-1'};
            app.DatasetDropDown.BackgroundColor = [1 1 1];
            app.DatasetDropDown.Position = [71 336 129 22];
            app.DatasetDropDown.Value = 'SENTINEL-1';

            % Create StartdateDatePickerLabel
            app.StartdateDatePickerLabel = uilabel(app.FilterPanel);
            app.StartdateDatePickerLabel.HorizontalAlignment = 'right';
            app.StartdateDatePickerLabel.FontColor = [0.149 0.149 0.149];
            app.StartdateDatePickerLabel.Position = [368 186 57 22];
            app.StartdateDatePickerLabel.Text = 'Start date';

            % Create StartdateDatePicker
            app.StartdateDatePicker = uidatepicker(app.FilterPanel);
            app.StartdateDatePicker.Limits = [datetime([0 1 1]) datetime([2026 12 31])];
            app.StartdateDatePicker.DisplayFormat = 'dd/MMM/uuuu';
            app.StartdateDatePicker.FontColor = [0.149 0.149 0.149];
            app.StartdateDatePicker.Placeholder = 'dd/mm/yyyy';
            app.StartdateDatePicker.Position = [440 186 104 22];

            % Create EnddateDatePickerLabel
            app.EnddateDatePickerLabel = uilabel(app.FilterPanel);
            app.EnddateDatePickerLabel.HorizontalAlignment = 'right';
            app.EnddateDatePickerLabel.FontColor = [0.149 0.149 0.149];
            app.EnddateDatePickerLabel.Position = [607 187 53 22];
            app.EnddateDatePickerLabel.Text = 'End date';

            % Create EnddateDatePicker
            app.EnddateDatePicker = uidatepicker(app.FilterPanel);
            app.EnddateDatePicker.Limits = [datetime([0 1 1]) datetime([2026 12 31])];
            app.EnddateDatePicker.DisplayFormat = 'dd/MMM/uuuu';
            app.EnddateDatePicker.FontColor = [0.149 0.149 0.149];
            app.EnddateDatePicker.Placeholder = 'dd/mm/yyyy';
            app.EnddateDatePicker.Position = [675 187 102 22];

            % Create ClearPathFrameButton
            app.ClearPathFrameButton = uibutton(app.FilterPanel, 'push');
            app.ClearPathFrameButton.ButtonPushedFcn = createCallbackFcn(app, @ClearPathFrameButtonPushed, true);
            app.ClearPathFrameButton.BackgroundColor = [1 1 1];
            app.ClearPathFrameButton.FontColor = [0.149 0.149 0.149];
            app.ClearPathFrameButton.Position = [375 95 116 22];
            app.ClearPathFrameButton.Text = 'Clear Path / Frame';

            % Create FileTypeLabel
            app.FileTypeLabel = uilabel(app.FilterPanel);
            app.FileTypeLabel.FontColor = [0.149 0.149 0.149];
            app.FileTypeLabel.Position = [14 313 53 22];
            app.FileTypeLabel.Text = 'File Type';

            % Create L1DetectedHighResDualPolGRDHDCheckBox
            app.L1DetectedHighResDualPolGRDHDCheckBox = uicheckbox(app.FilterPanel);
            app.L1DetectedHighResDualPolGRDHDCheckBox.Text = 'L1 Detected High-Res Dual-Pol (GRD-HD)';
            app.L1DetectedHighResDualPolGRDHDCheckBox.FontColor = [0.149 0.149 0.149];
            app.L1DetectedHighResDualPolGRDHDCheckBox.Position = [14 296 250 22];

            % Create L1DetectedMidResDualPolGRDMDCheckBox
            app.L1DetectedMidResDualPolGRDMDCheckBox = uicheckbox(app.FilterPanel);
            app.L1DetectedMidResDualPolGRDMDCheckBox.Text = 'L1 Detected Mid-Res Dual-Pol (GRD-MD)';
            app.L1DetectedMidResDualPolGRDMDCheckBox.FontColor = [0.149 0.149 0.149];
            app.L1DetectedMidResDualPolGRDMDCheckBox.Position = [14 279 246 22];

            % Create L1DetectedMidResSinglePolGRDMSCheckBox
            app.L1DetectedMidResSinglePolGRDMSCheckBox = uicheckbox(app.FilterPanel);
            app.L1DetectedMidResSinglePolGRDMSCheckBox.Text = 'L1 Detected Mid-Res Single-Pol (GRD-MS)';
            app.L1DetectedMidResSinglePolGRDMSCheckBox.FontColor = [0.149 0.149 0.149];
            app.L1DetectedMidResSinglePolGRDMSCheckBox.Position = [14 263 254 22];

            % Create L1DetectedHighResSinglePolGRDHSCheckBox
            app.L1DetectedHighResSinglePolGRDHSCheckBox = uicheckbox(app.FilterPanel);
            app.L1DetectedHighResSinglePolGRDHSCheckBox.Text = 'L1 Detected High-Res Single-Pol (GRD-HS)';
            app.L1DetectedHighResSinglePolGRDHSCheckBox.FontColor = [0.149 0.149 0.149];
            app.L1DetectedHighResSinglePolGRDHSCheckBox.Position = [14 245 259 22];

            % Create L1SingleLookComplexSLCCheckBox
            app.L1SingleLookComplexSLCCheckBox = uicheckbox(app.FilterPanel);
            app.L1SingleLookComplexSLCCheckBox.Text = 'L1 Single Look Complex (SLC)';
            app.L1SingleLookComplexSLCCheckBox.FontColor = [0.149 0.149 0.149];
            app.L1SingleLookComplexSLCCheckBox.Position = [14 228 186 22];
            app.L1SingleLookComplexSLCCheckBox.Value = true;

            % Create L2OceanOCNCheckBox
            app.L2OceanOCNCheckBox = uicheckbox(app.FilterPanel);
            app.L2OceanOCNCheckBox.Text = 'L2 Ocean (OCN)';
            app.L2OceanOCNCheckBox.FontColor = [0.149 0.149 0.149];
            app.L2OceanOCNCheckBox.Position = [14 211 111 22];

            % Create L0RawDataRAWCheckBox
            app.L0RawDataRAWCheckBox = uicheckbox(app.FilterPanel);
            app.L0RawDataRAWCheckBox.Text = 'L0 Raw Data (RAW)';
            app.L0RawDataRAWCheckBox.FontColor = [0.149 0.149 0.149];
            app.L0RawDataRAWCheckBox.Position = [14 194 129 22];

            % Create XMLMetadataGRDMSCheckBox
            app.XMLMetadataGRDMSCheckBox = uicheckbox(app.FilterPanel);
            app.XMLMetadataGRDMSCheckBox.Text = 'XML Metadata (GRD-MS)';
            app.XMLMetadataGRDMSCheckBox.FontColor = [0.149 0.149 0.149];
            app.XMLMetadataGRDMSCheckBox.Position = [14 177 159 22];

            % Create XMLMetadataGRDHDCheckBox
            app.XMLMetadataGRDHDCheckBox = uicheckbox(app.FilterPanel);
            app.XMLMetadataGRDHDCheckBox.Text = 'XML Metadata (GRD-HD)';
            app.XMLMetadataGRDHDCheckBox.FontColor = [0.149 0.149 0.149];
            app.XMLMetadataGRDHDCheckBox.Position = [14 160 158 22];

            % Create XMLMetadataRAWCheckBox
            app.XMLMetadataRAWCheckBox = uicheckbox(app.FilterPanel);
            app.XMLMetadataRAWCheckBox.Text = 'XML Metadata (RAW)';
            app.XMLMetadataRAWCheckBox.FontColor = [0.149 0.149 0.149];
            app.XMLMetadataRAWCheckBox.Position = [14 143 138 22];

            % Create XMLMetadataGRDHSCheckBox
            app.XMLMetadataGRDHSCheckBox = uicheckbox(app.FilterPanel);
            app.XMLMetadataGRDHSCheckBox.Text = 'XML Metadata (GRD-HS)';
            app.XMLMetadataGRDHSCheckBox.FontColor = [0.149 0.149 0.149];
            app.XMLMetadataGRDHSCheckBox.Position = [14 126 157 22];

            % Create XMLMetadataGRDMDCheckBox
            app.XMLMetadataGRDMDCheckBox = uicheckbox(app.FilterPanel);
            app.XMLMetadataGRDMDCheckBox.Text = 'XML Metadata (GRD-MD)';
            app.XMLMetadataGRDMDCheckBox.FontColor = [0.149 0.149 0.149];
            app.XMLMetadataGRDMDCheckBox.Position = [14 109 159 22];

            % Create XMLMetadataSLCCheckBox
            app.XMLMetadataSLCCheckBox = uicheckbox(app.FilterPanel);
            app.XMLMetadataSLCCheckBox.Text = 'XML Metadata (SLC)';
            app.XMLMetadataSLCCheckBox.FontColor = [0.149 0.149 0.149];
            app.XMLMetadataSLCCheckBox.Position = [14 92 133 22];

            % Create XMLMetadataOCNCheckBox
            app.XMLMetadataOCNCheckBox = uicheckbox(app.FilterPanel);
            app.XMLMetadataOCNCheckBox.Text = 'XML Metadata (OCN)';
            app.XMLMetadataOCNCheckBox.FontColor = [0.149 0.149 0.149];
            app.XMLMetadataOCNCheckBox.Position = [14 75 137 22];

            % Create OrbitDirectionLabel
            app.OrbitDirectionLabel = uilabel(app.FilterPanel);
            app.OrbitDirectionLabel.FontColor = [0.149 0.149 0.149];
            app.OrbitDirectionLabel.Position = [14 54 82 22];
            app.OrbitDirectionLabel.Text = 'Orbit Direction';

            % Create ASCENDINGCheckBox
            app.ASCENDINGCheckBox = uicheckbox(app.FilterPanel);
            app.ASCENDINGCheckBox.Text = 'ASCENDING';
            app.ASCENDINGCheckBox.FontColor = [0.149 0.149 0.149];
            app.ASCENDINGCheckBox.Position = [14 37 92 22];

            % Create DESCENDINGCheckBox
            app.DESCENDINGCheckBox = uicheckbox(app.FilterPanel);
            app.DESCENDINGCheckBox.Text = 'DESCENDING';
            app.DESCENDINGCheckBox.FontColor = [0.149 0.149 0.149];
            app.DESCENDINGCheckBox.Position = [14 20 101 22];

            % Create BeamModeLabel
            app.BeamModeLabel = uilabel(app.FilterPanel);
            app.BeamModeLabel.FontColor = [0.149 0.149 0.149];
            app.BeamModeLabel.Position = [288 313 70 22];
            app.BeamModeLabel.Text = 'Beam Mode';

            % Create IWCheckBox
            app.IWCheckBox = uicheckbox(app.FilterPanel);
            app.IWCheckBox.Text = 'IW';
            app.IWCheckBox.FontColor = [0.149 0.149 0.149];
            app.IWCheckBox.Position = [288 296 36 22];
            app.IWCheckBox.Value = true;

            % Create EWCheckBox
            app.EWCheckBox = uicheckbox(app.FilterPanel);
            app.EWCheckBox.Text = 'EW';
            app.EWCheckBox.FontColor = [0.149 0.149 0.149];
            app.EWCheckBox.Position = [288 279 40 22];

            % Create S1CheckBox
            app.S1CheckBox = uicheckbox(app.FilterPanel);
            app.S1CheckBox.Text = 'S1';
            app.S1CheckBox.FontColor = [0.149 0.149 0.149];
            app.S1CheckBox.Position = [288 262 36 22];

            % Create S2CheckBox
            app.S2CheckBox = uicheckbox(app.FilterPanel);
            app.S2CheckBox.Text = 'S2';
            app.S2CheckBox.FontColor = [0.149 0.149 0.149];
            app.S2CheckBox.Position = [288 245 36 22];

            % Create S3CheckBox
            app.S3CheckBox = uicheckbox(app.FilterPanel);
            app.S3CheckBox.Text = 'S3';
            app.S3CheckBox.FontColor = [0.149 0.149 0.149];
            app.S3CheckBox.Position = [288 228 36 22];

            % Create S4CheckBox
            app.S4CheckBox = uicheckbox(app.FilterPanel);
            app.S4CheckBox.Text = 'S4';
            app.S4CheckBox.FontColor = [0.149 0.149 0.149];
            app.S4CheckBox.Position = [288 211 36 22];

            % Create S5CheckBox
            app.S5CheckBox = uicheckbox(app.FilterPanel);
            app.S5CheckBox.Text = 'S5';
            app.S5CheckBox.FontColor = [0.149 0.149 0.149];
            app.S5CheckBox.Position = [288 194 36 22];

            % Create S6CheckBox
            app.S6CheckBox = uicheckbox(app.FilterPanel);
            app.S6CheckBox.Text = 'S6';
            app.S6CheckBox.FontColor = [0.149 0.149 0.149];
            app.S6CheckBox.Position = [288 177 36 22];

            % Create WVCheckBox
            app.WVCheckBox = uicheckbox(app.FilterPanel);
            app.WVCheckBox.Text = 'WV';
            app.WVCheckBox.FontColor = [0.149 0.149 0.149];
            app.WVCheckBox.Position = [288 160 40 22];

            % Create PolarizationLabel
            app.PolarizationLabel = uilabel(app.FilterPanel);
            app.PolarizationLabel.FontColor = [0.149 0.149 0.149];
            app.PolarizationLabel.Position = [288 139 68 22];
            app.PolarizationLabel.Text = 'Polarization';

            % Create VVCheckBox
            app.VVCheckBox = uicheckbox(app.FilterPanel);
            app.VVCheckBox.Text = 'VV';
            app.VVCheckBox.FontColor = [0.149 0.149 0.149];
            app.VVCheckBox.Position = [288 121 37 22];

            % Create HHCheckBox
            app.HHCheckBox = uicheckbox(app.FilterPanel);
            app.HHCheckBox.Text = 'HH';
            app.HHCheckBox.FontColor = [0.149 0.149 0.149];
            app.HHCheckBox.Position = [288 104 38 22];

            % Create VVHHCheckBox
            app.VVHHCheckBox = uicheckbox(app.FilterPanel);
            app.VVHHCheckBox.Text = 'VV+HH';
            app.VVHHCheckBox.FontColor = [0.149 0.149 0.149];
            app.VVHHCheckBox.Position = [288 87 61 22];

            % Create HHHVCheckBox
            app.HHHVCheckBox = uicheckbox(app.FilterPanel);
            app.HHHVCheckBox.Text = 'HH+HV';
            app.HHHVCheckBox.FontColor = [0.149 0.149 0.149];
            app.HHHVCheckBox.Position = [288 70 62 22];

            % Create DualHHCheckBox
            app.DualHHCheckBox = uicheckbox(app.FilterPanel);
            app.DualHHCheckBox.Text = 'Dual HH';
            app.DualHHCheckBox.FontColor = [0.149 0.149 0.149];
            app.DualHHCheckBox.Position = [288 53 66 22];

            % Create DualHVCheckBox
            app.DualHVCheckBox = uicheckbox(app.FilterPanel);
            app.DualHVCheckBox.Text = 'Dual HV';
            app.DualHVCheckBox.FontColor = [0.149 0.149 0.149];
            app.DualHVCheckBox.Position = [288 36 66 22];

            % Create DualVHCheckBox
            app.DualVHCheckBox = uicheckbox(app.FilterPanel);
            app.DualVHCheckBox.Text = 'Dual VH';
            app.DualVHCheckBox.FontColor = [0.149 0.149 0.149];
            app.DualVHCheckBox.Position = [288 19 66 22];

            % Create DualVVCheckBox
            app.DualVVCheckBox = uicheckbox(app.FilterPanel);
            app.DualVVCheckBox.Text = 'Dual VV';
            app.DualVVCheckBox.FontColor = [0.149 0.149 0.149];
            app.DualVVCheckBox.Position = [288 2 65 22];

            % Create SubtypeLabel
            app.SubtypeLabel = uilabel(app.FilterPanel);
            app.SubtypeLabel.FontColor = [0.149 0.149 0.149];
            app.SubtypeLabel.Position = [376 313 49 22];
            app.SubtypeLabel.Text = 'Subtype';

            % Create SACheckBox
            app.SACheckBox = uicheckbox(app.FilterPanel);
            app.SACheckBox.Text = 'SA';
            app.SACheckBox.FontColor = [0.149 0.149 0.149];
            app.SACheckBox.Position = [375 296 37 22];

            % Create SBCheckBox
            app.SBCheckBox = uicheckbox(app.FilterPanel);
            app.SBCheckBox.Text = 'SB';
            app.SBCheckBox.FontColor = [0.149 0.149 0.149];
            app.SBCheckBox.Position = [375 279 37 22];

            % Create SCCheckBox
            app.SCCheckBox = uicheckbox(app.FilterPanel);
            app.SCCheckBox.Text = 'SC';
            app.SCCheckBox.FontColor = [0.149 0.149 0.149];
            app.SCCheckBox.Position = [375 262 38 22];

            % Create SDCheckBox
            app.SDCheckBox = uicheckbox(app.FilterPanel);
            app.SDCheckBox.Text = 'SD';
            app.SDCheckBox.FontColor = [0.149 0.149 0.149];
            app.SDCheckBox.Position = [375 245 38 22];

            % Create GroupIDEditFieldLabel
            app.GroupIDEditFieldLabel = uilabel(app.FilterPanel);
            app.GroupIDEditFieldLabel.HorizontalAlignment = 'right';
            app.GroupIDEditFieldLabel.FontColor = [0.149 0.149 0.149];
            app.GroupIDEditFieldLabel.Position = [443 313 54 22];
            app.GroupIDEditFieldLabel.Text = 'Group ID';

            % Create GroupIDEditField
            app.GroupIDEditField = uieditfield(app.FilterPanel, 'numeric');
            app.GroupIDEditField.AllowEmpty = 'on';
            app.GroupIDEditField.FontColor = [0.149 0.149 0.149];
            app.GroupIDEditField.Position = [512 313 100 22];
            app.GroupIDEditField.Value = [];

            % Create ResetallfiltersButton
            app.ResetallfiltersButton = uibutton(app.FilterPanel, 'push');
            app.ResetallfiltersButton.ButtonPushedFcn = createCallbackFcn(app, @ResetallfiltersButtonPushed, true);
            app.ResetallfiltersButton.Position = [652 313 100 22];
            app.ResetallfiltersButton.Text = 'Reset all filters';

            % Create SamplingRateEditFieldLabel
            app.SamplingRateEditFieldLabel = uilabel(app.FilterPanel);
            app.SamplingRateEditFieldLabel.HorizontalAlignment = 'right';
            app.SamplingRateEditFieldLabel.FontColor = [0.149 0.149 0.149];
            app.SamplingRateEditFieldLabel.Position = [443 275 84 22];
            app.SamplingRateEditFieldLabel.Text = 'Sampling Rate';

            % Create SamplingRateEditField
            app.SamplingRateEditField = uieditfield(app.FilterPanel, 'numeric');
            app.SamplingRateEditField.FontColor = [0.149 0.149 0.149];
            app.SamplingRateEditField.Position = [542 275 71 22];

            % Create SamplingUnitDropDown
            app.SamplingUnitDropDown = uidropdown(app.FilterPanel);
            app.SamplingUnitDropDown.Items = {'Week', 'Month', 'Year'};
            app.SamplingUnitDropDown.FontColor = [0 0 0];
            app.SamplingUnitDropDown.BackgroundColor = [1 1 1];
            app.SamplingUnitDropDown.Position = [622 275 100 22];
            app.SamplingUnitDropDown.Value = 'Month';

            % Create MinLongitudeEditFieldLabel
            app.MinLongitudeEditFieldLabel = uilabel(app.DownloadTab);
            app.MinLongitudeEditFieldLabel.HorizontalAlignment = 'right';
            app.MinLongitudeEditFieldLabel.Position = [964 111 80 22];
            app.MinLongitudeEditFieldLabel.Text = 'Min Longitude';

            % Create MinLongitudeEditField
            app.MinLongitudeEditField = uieditfield(app.DownloadTab, 'numeric');
            app.MinLongitudeEditField.ValueChangedFcn = createCallbackFcn(app, @MinLongitudeEditFieldValueChanged, true);
            app.MinLongitudeEditField.Position = [1059 111 100 22];

            % Create MaxLongitudeEditFieldLabel
            app.MaxLongitudeEditFieldLabel = uilabel(app.DownloadTab);
            app.MaxLongitudeEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxLongitudeEditFieldLabel.Position = [960 85 84 22];
            app.MaxLongitudeEditFieldLabel.Text = 'Max Longitude';

            % Create MaxLongitudeEditField
            app.MaxLongitudeEditField = uieditfield(app.DownloadTab, 'numeric');
            app.MaxLongitudeEditField.ValueChangedFcn = createCallbackFcn(app, @MaxLongitudeEditFieldValueChanged, true);
            app.MaxLongitudeEditField.Position = [1059 85 100 22];

            % Create MinLatitudeEditFieldLabel
            app.MinLatitudeEditFieldLabel = uilabel(app.DownloadTab);
            app.MinLatitudeEditFieldLabel.HorizontalAlignment = 'right';
            app.MinLatitudeEditFieldLabel.Position = [974 59 70 22];
            app.MinLatitudeEditFieldLabel.Text = 'Min Latitude';

            % Create MinLatitudeEditField
            app.MinLatitudeEditField = uieditfield(app.DownloadTab, 'numeric');
            app.MinLatitudeEditField.ValueChangedFcn = createCallbackFcn(app, @MinLatitudeEditFieldValueChanged, true);
            app.MinLatitudeEditField.Position = [1059 59 100 22];

            % Create MaxLatitudeEditFieldLabel
            app.MaxLatitudeEditFieldLabel = uilabel(app.DownloadTab);
            app.MaxLatitudeEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxLatitudeEditFieldLabel.Position = [970 34 74 22];
            app.MaxLatitudeEditFieldLabel.Text = 'Max Latitude';

            % Create MaxLatitudeEditField
            app.MaxLatitudeEditField = uieditfield(app.DownloadTab, 'numeric');
            app.MaxLatitudeEditField.ValueChangedFcn = createCallbackFcn(app, @MaxLatitudeEditFieldValueChanged, true);
            app.MaxLatitudeEditField.Position = [1059 34 100 22];

            % Create EarthdataUsernameEditField
            app.EarthdataUsernameEditField = uieditfield(app.DownloadTab, 'text');
            app.EarthdataUsernameEditField.Placeholder = 'Username';
            app.EarthdataUsernameEditField.Position = [949 455 100 22];

            % Create EarthdataPasswordEditField
            app.EarthdataPasswordEditField = uieditfield(app.DownloadTab, 'text');
            app.EarthdataPasswordEditField.Placeholder = 'Password';
            app.EarthdataPasswordEditField.Position = [949 423 100 22];

            % Create LoginButton
            app.LoginButton = uibutton(app.DownloadTab, 'push');
            app.LoginButton.ButtonPushedFcn = createCallbackFcn(app, @LoginButtonPushed, true);
            app.LoginButton.Position = [1059 422 66 22];
            app.LoginButton.Text = 'Login';

            % Create LoginFeedbackLabel
            app.LoginFeedbackLabel = uilabel(app.DownloadTab);
            app.LoginFeedbackLabel.Position = [949 391 239 22];
            app.LoginFeedbackLabel.Text = '';

            % Create SignedInLabel
            app.SignedInLabel = uilabel(app.DownloadTab);
            app.SignedInLabel.Position = [946 455 239 22];
            app.SignedInLabel.Text = '';

            % Create SignOutButton
            app.SignOutButton = uibutton(app.DownloadTab, 'push');
            app.SignOutButton.ButtonPushedFcn = createCallbackFcn(app, @SignOutButtonPushed, true);
            app.SignOutButton.Position = [1124 455 58 22];
            app.SignOutButton.Text = 'Sign Out';

            % Create DownloaderToolbarPanel
            app.DownloaderToolbarPanel = uipanel(app.DownloadTab);
            app.DownloaderToolbarPanel.Position = [3 448 912 30];

            % Create ShowFiltersButton
            app.ShowFiltersButton = uibutton(app.DownloaderToolbarPanel, 'push');
            app.ShowFiltersButton.ButtonPushedFcn = createCallbackFcn(app, @ShowFiltersButtonPushed, true);
            app.ShowFiltersButton.Position = [3 0 100 29];
            app.ShowFiltersButton.Text = 'Show Filters';

            % Create DrawRectangleButton
            app.DrawRectangleButton = uibutton(app.DownloaderToolbarPanel, 'push');
            app.DrawRectangleButton.ButtonPushedFcn = createCallbackFcn(app, @DrawRectangleButtonPushed, true);
            app.DrawRectangleButton.Position = [102 0 100 29];
            app.DrawRectangleButton.Text = 'Draw Rectangle';

            % Create DrawPolygonButton
            app.DrawPolygonButton = uibutton(app.DownloaderToolbarPanel, 'push');
            app.DrawPolygonButton.ButtonPushedFcn = createCallbackFcn(app, @DrawPolygonButtonPushed, true);
            app.DrawPolygonButton.Position = [202 0 100 29];
            app.DrawPolygonButton.Text = 'Draw Polygon';

            % Create LoadShapefileButton
            app.LoadShapefileButton = uibutton(app.DownloaderToolbarPanel, 'push');
            app.LoadShapefileButton.Position = [301 0 100 29];
            app.LoadShapefileButton.Text = 'Load Shapefile';

            % Create ClearAOIButton
            app.ClearAOIButton = uibutton(app.DownloaderToolbarPanel, 'push');
            app.ClearAOIButton.ButtonPushedFcn = createCallbackFcn(app, @ClearAOIButtonPushed, true);
            app.ClearAOIButton.Position = [408 0 100 29];
            app.ClearAOIButton.Text = 'Clear AOI';

            % Create SearchASFButton
            app.SearchASFButton = uibutton(app.DownloaderToolbarPanel, 'push');
            app.SearchASFButton.ButtonPushedFcn = createCallbackFcn(app, @SearchASFButtonPushed, true);
            app.SearchASFButton.Position = [746 1 100 29];
            app.SearchASFButton.Text = 'Search ASF';

            % Create LoadLastDownloadButton
            app.LoadLastDownloadButton = uibutton(app.DownloaderToolbarPanel, 'push');
            app.LoadLastDownloadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadLastDownloadButtonPushed, true);
            app.LoadLastDownloadButton.Position = [520 0 124 29];
            app.LoadLastDownloadButton.Text = 'Load Last Download';

            % Create TogglePreviewButton
            app.TogglePreviewButton = uibutton(app.DownloaderToolbarPanel, 'push');
            app.TogglePreviewButton.ButtonPushedFcn = createCallbackFcn(app, @TogglePreviewButtonPushed, true);
            app.TogglePreviewButton.Position = [644 0 100 29];
            app.TogglePreviewButton.Text = 'Show Preview';

            % Create PreviewPanel
            app.PreviewPanel = uipanel(app.DownloadTab);
            app.PreviewPanel.Title = 'Preview';
            app.PreviewPanel.BackgroundColor = [1 1 1];
            app.PreviewPanel.Position = [809 23 387 421];

            % Create PreviewResultsTable
            app.PreviewResultsTable = uitable(app.PreviewPanel);
            app.PreviewResultsTable.ColumnName = {'Select'; 'Date'; 'Time'; 'Path'; 'Frame'; 'Direction'; 'Size'};
            app.PreviewResultsTable.RowName = {};
            app.PreviewResultsTable.CellEditCallback = createCallbackFcn(app, @PreviewResultsTableCellSelection, true);
            app.PreviewResultsTable.Position = [0 187 380 185];

            % Create DownloadSelectedButton
            app.DownloadSelectedButton = uibutton(app.PreviewPanel, 'push');
            app.DownloadSelectedButton.ButtonPushedFcn = createCallbackFcn(app, @DownloadSelectedButtonPushed, true);
            app.DownloadSelectedButton.Position = [37 22 118 22];
            app.DownloadSelectedButton.Text = 'Download Selected';

            % Create PreviewSelectedLabel
            app.PreviewSelectedLabel = uilabel(app.PreviewPanel);
            app.PreviewSelectedLabel.Position = [151 376 98 22];
            app.PreviewSelectedLabel.Text = 'Preview Selected';

            % Create SelectAllCheckBox
            app.SelectAllCheckBox = uicheckbox(app.PreviewPanel);
            app.SelectAllCheckBox.ValueChangedFcn = createCallbackFcn(app, @SelectAllCheckBoxValueChanged, true);
            app.SelectAllCheckBox.Text = 'Select All';
            app.SelectAllCheckBox.Position = [29 376 71 22];

            % Create DownloadAllButton
            app.DownloadAllButton = uibutton(app.PreviewPanel, 'push');
            app.DownloadAllButton.ButtonPushedFcn = createCallbackFcn(app, @DownloadAllButtonPushed, true);
            app.DownloadAllButton.Position = [169 22 100 22];
            app.DownloadAllButton.Text = 'Download All';

            % Create ToggleFootprintsButton
            app.ToggleFootprintsButton = uibutton(app.PreviewPanel, 'push');
            app.ToggleFootprintsButton.ButtonPushedFcn = createCallbackFcn(app, @ToggleFootprintsButtonPushed, true);
            app.ToggleFootprintsButton.Position = [15 157 101 22];
            app.ToggleFootprintsButton.Text = 'Show Footprints';

            % Create RecommendedFramePathLabel
            app.RecommendedFramePathLabel = uilabel(app.PreviewPanel);
            app.RecommendedFramePathLabel.Position = [23 125 301 22];
            app.RecommendedFramePathLabel.Text = '';

            % Create DownloadProgressGauge
            app.DownloadProgressGauge = uigauge(app.PreviewPanel, 'linear');
            app.DownloadProgressGauge.MinorTicks = [0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80 82 84 86 88 90 92 94 96 98 100];
            app.DownloadProgressGauge.Position = [23 50 203 41];

            % Create DownloadProgressLabel
            app.DownloadProgressLabel = uilabel(app.PreviewPanel);
            app.DownloadProgressLabel.Position = [23 98 344 22];
            app.DownloadProgressLabel.Text = '';

            % Create GlobalVariablesTab_3
            app.GlobalVariablesTab_3 = uitab(app.TabGroup);
            app.GlobalVariablesTab_3.Title = 'Global Variables';
            app.GlobalVariablesTab_3.BackgroundColor = [1 1 1];

            % Create PythonEnvironmentLabel
            app.PythonEnvironmentLabel = uilabel(app.GlobalVariablesTab_3);
            app.PythonEnvironmentLabel.HorizontalAlignment = 'right';
            app.PythonEnvironmentLabel.Position = [43 358 161 22];
            app.PythonEnvironmentLabel.Text = 'Custom Python Environment:';

            % Create UpdateImagesTable
            app.UpdateImagesTable = uitable(app.GlobalVariablesTab_3);
            app.UpdateImagesTable.ColumnName = {'Use'; 'Scene'; 'Acquisition'; 'Platform'; 'Polarization'};
            app.UpdateImagesTable.RowName = {};
            app.UpdateImagesTable.ColumnEditable = [true false false false false];
            app.UpdateImagesTable.Visible = 'off';
            app.UpdateImagesTable.Position = [520 90 520 125];

            % Create DeselectAllUpdateImagesButton
            app.DeselectAllUpdateImagesButton = uibutton(app.GlobalVariablesTab_3, 'push');
            app.DeselectAllUpdateImagesButton.ButtonPushedFcn = createCallbackFcn(app, @DeselectAllUpdateImagesButtonPushed, true);
            app.DeselectAllUpdateImagesButton.BackgroundColor = [1 1 1];
            app.DeselectAllUpdateImagesButton.Visible = 'off';
            app.DeselectAllUpdateImagesButton.Position = [866 221 90 23];
            app.DeselectAllUpdateImagesButton.Text = 'Deselect all';

            % Create SelectAllUpdateImagesButton
            app.SelectAllUpdateImagesButton = uibutton(app.GlobalVariablesTab_3, 'push');
            app.SelectAllUpdateImagesButton.ButtonPushedFcn = createCallbackFcn(app, @SelectAllUpdateImagesButtonPushed, true);
            app.SelectAllUpdateImagesButton.BackgroundColor = [1 1 1];
            app.SelectAllUpdateImagesButton.Visible = 'off';
            app.SelectAllUpdateImagesButton.Position = [775 221 80 23];
            app.SelectAllUpdateImagesButton.Text = 'Select all';

            % Create UpdateImagesTableLabel
            app.UpdateImagesTableLabel = uilabel(app.GlobalVariablesTab_3);
            app.UpdateImagesTableLabel.FontWeight = 'bold';
            app.UpdateImagesTableLabel.Visible = 'off';
            app.UpdateImagesTableLabel.Position = [520 221 92 22];
            app.UpdateImagesTableLabel.Text = 'Images found';

            % Create UpdateSearchParametersTextArea
            app.UpdateSearchParametersTextArea = uitextarea(app.GlobalVariablesTab_3);
            app.UpdateSearchParametersTextArea.Editable = 'off';
            app.UpdateSearchParametersTextArea.Visible = 'off';
            app.UpdateSearchParametersTextArea.Position = [64 90 420 125];
            app.UpdateSearchParametersTextArea.Value = {'File type: L1 Single Look Complex (SLC)'};

            % Create UpdateSearchParametersLabel
            app.UpdateSearchParametersLabel = uilabel(app.GlobalVariablesTab_3);
            app.UpdateSearchParametersLabel.FontWeight = 'bold';
            app.UpdateSearchParametersLabel.Visible = 'off';
            app.UpdateSearchParametersLabel.Position = [64 221 130 22];
            app.UpdateSearchParametersLabel.Text = 'Search parameters';

            % Create UpdateSearchStatusLabel
            app.UpdateSearchStatusLabel = uilabel(app.GlobalVariablesTab_3);
            app.UpdateSearchStatusLabel.Visible = 'off';
            app.UpdateSearchStatusLabel.Position = [64 249 700 22];
            app.UpdateSearchStatusLabel.Text = 'Click Search Images to query ASF.';

            % Create DownloadRunUpdateImagesButton
            app.DownloadRunUpdateImagesButton = uibutton(app.GlobalVariablesTab_3, 'push');
            app.DownloadRunUpdateImagesButton.ButtonPushedFcn = createCallbackFcn(app, @DownloadRunUpdateImagesButtonPushed, true);
            app.DownloadRunUpdateImagesButton.BackgroundColor = [1 1 1];
            app.DownloadRunUpdateImagesButton.Enable = 'off';
            app.DownloadRunUpdateImagesButton.Visible = 'off';
            app.DownloadRunUpdateImagesButton.Position = [755 278 178 24];
            app.DownloadRunUpdateImagesButton.Text = 'Download selected and run';

            % Create DownloadUpdateImagesButton
            app.DownloadUpdateImagesButton = uibutton(app.GlobalVariablesTab_3, 'push');
            app.DownloadUpdateImagesButton.ButtonPushedFcn = createCallbackFcn(app, @DownloadUpdateImagesButtonPushed, true);
            app.DownloadUpdateImagesButton.BackgroundColor = [1 1 1];
            app.DownloadUpdateImagesButton.Enable = 'off';
            app.DownloadUpdateImagesButton.Visible = 'off';
            app.DownloadUpdateImagesButton.Position = [610 278 132 24];
            app.DownloadUpdateImagesButton.Text = 'Download Selected';

            % Create SearchUpdateImagesButton
            app.SearchUpdateImagesButton = uibutton(app.GlobalVariablesTab_3, 'push');
            app.SearchUpdateImagesButton.ButtonPushedFcn = createCallbackFcn(app, @SearchUpdateImagesButtonPushed, true);
            app.SearchUpdateImagesButton.BackgroundColor = [1 1 1];
            app.SearchUpdateImagesButton.Enable = 'off';
            app.SearchUpdateImagesButton.Visible = 'off';
            app.SearchUpdateImagesButton.Position = [480 278 115 24];
            app.SearchUpdateImagesButton.Text = 'Search Images';

            % Create UpdateSearchEndDateDatePicker
            app.UpdateSearchEndDateDatePicker = uidatepicker(app.GlobalVariablesTab_3);
            app.UpdateSearchEndDateDatePicker.ValueChangedFcn = createCallbackFcn(app, @UpdateSearchEndDateDatePickerValueChanged, true);
            app.UpdateSearchEndDateDatePicker.Visible = 'off';
            app.UpdateSearchEndDateDatePicker.Position = [314 279 150 22];

            % Create UpdateSearchFromLabel
            app.UpdateSearchFromLabel = uilabel(app.GlobalVariablesTab_3);
            app.UpdateSearchFromLabel.HorizontalAlignment = 'right';
            app.UpdateSearchFromLabel.Visible = 'off';
            app.UpdateSearchFromLabel.Position = [44 279 255 22];
            app.UpdateSearchFromLabel.Text = 'Search new images from - to';

            % Create UpdateAlreadyProcessedDataCheckBox
            app.UpdateAlreadyProcessedDataCheckBox = uicheckbox(app.GlobalVariablesTab_3);
            app.UpdateAlreadyProcessedDataCheckBox.ValueChangedFcn = createCallbackFcn(app, @UpdateAlreadyProcessedDataCheckBoxValueChanged, true);
            app.UpdateAlreadyProcessedDataCheckBox.Text = 'Update already processed data';
            app.UpdateAlreadyProcessedDataCheckBox.Position = [46 315 211 22];

            % Create CustomPythonEnvironmentEditField
            app.CustomPythonEnvironmentEditField = uieditfield(app.GlobalVariablesTab_3, 'text');
            app.CustomPythonEnvironmentEditField.ValueChangedFcn = createCallbackFcn(app, @CustomPythonEnvironmentEditFieldValueChanged, true);
            app.CustomPythonEnvironmentEditField.FontName = 'Manrope';
            app.CustomPythonEnvironmentEditField.Position = [219 358 161 22];
            app.CustomPythonEnvironmentEditField.Value = 'python';

            % Create Label_2
            app.Label_2 = uilabel(app.GlobalVariablesTab_3);
            app.Label_2.HorizontalAlignment = 'center';
            app.Label_2.FontName = 'Manrope';
            app.Label_2.FontWeight = 'bold';
            app.Label_2.Position = [39 430 277 22];
            app.Label_2.Text = 'Name of the python 3.x environment in your os';

            % Create PythonEnvironmentDropDownLabel
            app.PythonEnvironmentDropDownLabel = uilabel(app.GlobalVariablesTab_3);
            app.PythonEnvironmentDropDownLabel.HorizontalAlignment = 'right';
            app.PythonEnvironmentDropDownLabel.Position = [45 398 113 22];
            app.PythonEnvironmentDropDownLabel.Text = 'Python Environment';

            % Create PythonEnvironmentDropDown
            app.PythonEnvironmentDropDown = uidropdown(app.GlobalVariablesTab_3);
            app.PythonEnvironmentDropDown.Items = {'python', 'python3', 'python3.11', 'Other'};
            app.PythonEnvironmentDropDown.ValueChangedFcn = createCallbackFcn(app, @PythonEnvironmentDropDownValueChanged, true);
            app.PythonEnvironmentDropDown.FontName = 'Manrope';
            app.PythonEnvironmentDropDown.BackgroundColor = [1 1 1];
            app.PythonEnvironmentDropDown.Position = [173 398 100 22];
            app.PythonEnvironmentDropDown.Value = 'python';

            % Create AOITab_2
            app.AOITab_2 = uitab(app.TabGroup);
            app.AOITab_2.Title = 'AOI';
            app.AOITab_2.BackgroundColor = [1 1 1];

            % Create MinlongitudeEditFieldLabel
            app.MinlongitudeEditFieldLabel = uilabel(app.AOITab_2);
            app.MinlongitudeEditFieldLabel.HorizontalAlignment = 'right';
            app.MinlongitudeEditFieldLabel.Position = [900 133 76 22];
            app.MinlongitudeEditFieldLabel.Text = 'Min longitude';

            % Create MinlongitudeEditField
            app.MinlongitudeEditField = uieditfield(app.AOITab_2, 'numeric');
            app.MinlongitudeEditField.Limits = [-180 180];
            app.MinlongitudeEditField.ValueDisplayFormat = '%.3f';
            app.MinlongitudeEditField.ValueChangedFcn = createCallbackFcn(app, @MinlongitudeEditFieldValueChanged, true);
            app.MinlongitudeEditField.Editable = 'off';
            app.MinlongitudeEditField.FontName = 'Manrope';
            app.MinlongitudeEditField.Position = [995 133 100 22];
            app.MinlongitudeEditField.Value = -180;

            % Create MaxlongitudeEditFieldLabel
            app.MaxlongitudeEditFieldLabel = uilabel(app.AOITab_2);
            app.MaxlongitudeEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxlongitudeEditFieldLabel.Position = [900 100 80 22];
            app.MaxlongitudeEditFieldLabel.Text = 'Max longitude';

            % Create MaxlongitudeEditField
            app.MaxlongitudeEditField = uieditfield(app.AOITab_2, 'numeric');
            app.MaxlongitudeEditField.Limits = [-180 180];
            app.MaxlongitudeEditField.ValueDisplayFormat = '%.3f';
            app.MaxlongitudeEditField.ValueChangedFcn = createCallbackFcn(app, @MaxlongitudeEditFieldValueChanged, true);
            app.MaxlongitudeEditField.Editable = 'off';
            app.MaxlongitudeEditField.FontName = 'Manrope';
            app.MaxlongitudeEditField.Position = [995 100 100 22];
            app.MaxlongitudeEditField.Value = 180;

            % Create MinlatitudeEditFieldLabel
            app.MinlatitudeEditFieldLabel = uilabel(app.AOITab_2);
            app.MinlatitudeEditFieldLabel.HorizontalAlignment = 'right';
            app.MinlatitudeEditFieldLabel.Position = [900 66 66 22];
            app.MinlatitudeEditFieldLabel.Text = 'Min latitude';

            % Create MinlatitudeEditField
            app.MinlatitudeEditField = uieditfield(app.AOITab_2, 'numeric');
            app.MinlatitudeEditField.Limits = [-90 90];
            app.MinlatitudeEditField.ValueDisplayFormat = '%.3f';
            app.MinlatitudeEditField.ValueChangedFcn = createCallbackFcn(app, @MinlatitudeEditFieldValueChanged, true);
            app.MinlatitudeEditField.Editable = 'off';
            app.MinlatitudeEditField.FontName = 'Manrope';
            app.MinlatitudeEditField.Position = [995 66 100 22];
            app.MinlatitudeEditField.Value = -90;

            % Create MaxlatitudeEditFieldLabel
            app.MaxlatitudeEditFieldLabel = uilabel(app.AOITab_2);
            app.MaxlatitudeEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxlatitudeEditFieldLabel.Position = [900 34 70 22];
            app.MaxlatitudeEditFieldLabel.Text = 'Max latitude';

            % Create MaxlatitudeEditField
            app.MaxlatitudeEditField = uieditfield(app.AOITab_2, 'numeric');
            app.MaxlatitudeEditField.Limits = [-90 90];
            app.MaxlatitudeEditField.ValueDisplayFormat = '%.3f';
            app.MaxlatitudeEditField.ValueChangedFcn = createCallbackFcn(app, @MaxlatitudeEditFieldValueChanged, true);
            app.MaxlatitudeEditField.Editable = 'off';
            app.MaxlatitudeEditField.FontName = 'Manrope';
            app.MaxlatitudeEditField.Position = [995 34 100 22];
            app.MaxlatitudeEditField.Value = 90;

            % Create AOIboxboundariesLabel_3
            app.AOIboxboundariesLabel_3 = uilabel(app.AOITab_2);
            app.AOIboxboundariesLabel_3.FontName = 'Manrope';
            app.AOIboxboundariesLabel_3.FontWeight = 'bold';
            app.AOIboxboundariesLabel_3.Position = [904 169 117 22];
            app.AOIboxboundariesLabel_3.Text = 'AOI box boundaries';

            % Create ResetSENMapButton
            app.ResetSENMapButton = uibutton(app.AOITab_2, 'push');
            app.ResetSENMapButton.ButtonPushedFcn = createCallbackFcn(app, @ResetSENMapButtonPushed, true);
            app.ResetSENMapButton.BackgroundColor = [1 1 1];
            app.ResetSENMapButton.FontName = 'Manrope';
            app.ResetSENMapButton.Position = [900 447 149 23];
            app.ResetSENMapButton.Text = 'Reset Map';

            % Create SENMapPanel
            app.SENMapPanel = uipanel(app.AOITab_2);
            app.SENMapPanel.BorderType = 'none';
            app.SENMapPanel.BackgroundColor = [1 1 1];
            app.SENMapPanel.FontName = 'Manrope';
            app.SENMapPanel.Position = [26 23 822 447];

            % Create MasterProcessingTab
            app.MasterProcessingTab = uitab(app.TabGroup);
            app.MasterProcessingTab.Title = 'Master Processing';
            app.MasterProcessingTab.BackgroundColor = [1 1 1];

            % Create MasterdateDatePickerLabel
            app.MasterdateDatePickerLabel = uilabel(app.MasterProcessingTab);
            app.MasterdateDatePickerLabel.HorizontalAlignment = 'right';
            app.MasterdateDatePickerLabel.FontWeight = 'bold';
            app.MasterdateDatePickerLabel.Position = [40 437 72 22];
            app.MasterdateDatePickerLabel.Text = 'Master date';

            % Create MasterdateDatePicker
            app.MasterdateDatePicker = uidatepicker(app.MasterProcessingTab);
            app.MasterdateDatePicker.ValueChangedFcn = createCallbackFcn(app, @MasterdateDatePickerValueChanged, true);
            app.MasterdateDatePicker.FontName = 'Manrope';
            app.MasterdateDatePicker.Enable = 'off';
            app.MasterdateDatePicker.Position = [127 437 150 22];
            app.MasterdateDatePicker.Value = datetime([2020 3 5]);

            % Create MasterprocessingCheckBox
            app.MasterprocessingCheckBox = uicheckbox(app.MasterProcessingTab);
            app.MasterprocessingCheckBox.ValueChangedFcn = createCallbackFcn(app, @MasterprocessingCheckBoxValueChanged, true);
            app.MasterprocessingCheckBox.Text = 'Master processing';
            app.MasterprocessingCheckBox.FontName = 'Manrope';
            app.MasterprocessingCheckBox.FontWeight = 'bold';
            app.MasterprocessingCheckBox.Position = [45 386 130 22];
            app.MasterprocessingCheckBox.Value = true;

            % Create PolarisationDropDownLabel
            app.PolarisationDropDownLabel = uilabel(app.MasterProcessingTab);
            app.PolarisationDropDownLabel.HorizontalAlignment = 'right';
            app.PolarisationDropDownLabel.FontWeight = 'bold';
            app.PolarisationDropDownLabel.Position = [45 316 74 22];
            app.PolarisationDropDownLabel.Text = 'Polarisation';

            % Create PolarisationDropDown
            app.PolarisationDropDown = uidropdown(app.MasterProcessingTab);
            app.PolarisationDropDown.Items = {'VV', 'VH', 'HH'};
            app.PolarisationDropDown.ValueChangedFcn = createCallbackFcn(app, @PolarisationDropDownValueChanged, true);
            app.PolarisationDropDown.FontName = 'Manrope';
            app.PolarisationDropDown.FontWeight = 'bold';
            app.PolarisationDropDown.BackgroundColor = [1 1 1];
            app.PolarisationDropDown.Position = [134 316 100 22];
            app.PolarisationDropDown.Value = 'VV';

            % Create AutoMasterCheckBox
            app.AutoMasterCheckBox = uicheckbox(app.MasterProcessingTab);
            app.AutoMasterCheckBox.ValueChangedFcn = createCallbackFcn(app, @AutoMasterCheckBoxValueChanged, true);
            app.AutoMasterCheckBox.Text = 'Auto-select master';
            app.AutoMasterCheckBox.FontName = 'Manrope';
            app.AutoMasterCheckBox.Position = [338 437 129 22];
            app.AutoMasterCheckBox.Value = true;

            % Create SlavesProcessingTab
            app.SlavesProcessingTab = uitab(app.TabGroup);
            app.SlavesProcessingTab.Title = 'Slaves Processing';
            app.SlavesProcessingTab.BackgroundColor = [1 1 1];

            % Create SlavesremovalafterprocessingCheckBox
            app.SlavesremovalafterprocessingCheckBox = uicheckbox(app.SlavesProcessingTab);
            app.SlavesremovalafterprocessingCheckBox.ValueChangedFcn = createCallbackFcn(app, @SlavesremovalafterprocessingCheckBoxValueChanged, true);
            app.SlavesremovalafterprocessingCheckBox.Text = 'Slaves removal after processing';
            app.SlavesremovalafterprocessingCheckBox.FontName = 'Manrope';
            app.SlavesremovalafterprocessingCheckBox.FontWeight = 'bold';
            app.SlavesremovalafterprocessingCheckBox.Position = [48 437 207 22];

            % Create DEMinterferogramDropDownLabel
            app.DEMinterferogramDropDownLabel = uilabel(app.SlavesProcessingTab);
            app.DEMinterferogramDropDownLabel.HorizontalAlignment = 'right';
            app.DEMinterferogramDropDownLabel.FontWeight = 'bold';
            app.DEMinterferogramDropDownLabel.Position = [40 312 113 22];
            app.DEMinterferogramDropDownLabel.Text = 'DEM interferogram';

            % Create DEMinterferogramDropDown
            app.DEMinterferogramDropDown = uidropdown(app.SlavesProcessingTab);
            app.DEMinterferogramDropDown.Items = {'SRTM 1Sec HGT', 'SRTM 3Sec', 'Copernicus 30m Global DEM', 'Copernicus 90m Global DEM', 'CDEM', 'GETASSE30', 'External DEM'};
            app.DEMinterferogramDropDown.ValueChangedFcn = createCallbackFcn(app, @DEMinterferogramDropDownValueChanged, true);
            app.DEMinterferogramDropDown.FontName = 'Manrope';
            app.DEMinterferogramDropDown.FontWeight = 'bold';
            app.DEMinterferogramDropDown.BackgroundColor = [1 1 1];
            app.DEMinterferogramDropDown.Position = [168 312 198 22];
            app.DEMinterferogramDropDown.Value = 'SRTM 1Sec HGT';

            % Create DEMcoregistrationDropDownLabel
            app.DEMcoregistrationDropDownLabel = uilabel(app.SlavesProcessingTab);
            app.DEMcoregistrationDropDownLabel.HorizontalAlignment = 'right';
            app.DEMcoregistrationDropDownLabel.FontWeight = 'bold';
            app.DEMcoregistrationDropDownLabel.Position = [41 231 115 22];
            app.DEMcoregistrationDropDownLabel.Text = 'DEM coregistration';

            % Create DEMcoregistrationDropDown
            app.DEMcoregistrationDropDown = uidropdown(app.SlavesProcessingTab);
            app.DEMcoregistrationDropDown.Items = {'SRTM 1Sec HGT', 'SRTM 3Sec', 'Copernicus 30m Global DEM', 'Copernicus 90m Global DEM', 'CDEM', 'GETASSE30', 'External DEM'};
            app.DEMcoregistrationDropDown.ValueChangedFcn = createCallbackFcn(app, @DEMcoregistrationDropDownValueChanged, true);
            app.DEMcoregistrationDropDown.FontName = 'Manrope';
            app.DEMcoregistrationDropDown.FontWeight = 'bold';
            app.DEMcoregistrationDropDown.BackgroundColor = [1 1 1];
            app.DEMcoregistrationDropDown.Position = [166 231 201 22];
            app.DEMcoregistrationDropDown.Value = 'SRTM 1Sec HGT';

            % Create DEMifgpathEditFieldLabel
            app.DEMifgpathEditFieldLabel = uilabel(app.SlavesProcessingTab);
            app.DEMifgpathEditFieldLabel.HorizontalAlignment = 'right';
            app.DEMifgpathEditFieldLabel.Position = [51 270 74 22];
            app.DEMifgpathEditFieldLabel.Text = 'DEM ifg path';

            % Create DEMifgpathEditField
            app.DEMifgpathEditField = uieditfield(app.SlavesProcessingTab, 'text');
            app.DEMifgpathEditField.ValueChangedFcn = createCallbackFcn(app, @DEMifgpathEditFieldValueChanged, true);
            app.DEMifgpathEditField.HorizontalAlignment = 'right';
            app.DEMifgpathEditField.FontName = 'Manrope';
            app.DEMifgpathEditField.Position = [140 270 503 22];

            % Create DEMcoregpathEditFieldLabel
            app.DEMcoregpathEditFieldLabel = uilabel(app.SlavesProcessingTab);
            app.DEMcoregpathEditFieldLabel.HorizontalAlignment = 'right';
            app.DEMcoregpathEditFieldLabel.Position = [51 190 92 22];
            app.DEMcoregpathEditFieldLabel.Text = 'DEM coreg path';

            % Create DEMcoregpathEditField
            app.DEMcoregpathEditField = uieditfield(app.SlavesProcessingTab, 'text');
            app.DEMcoregpathEditField.ValueChangedFcn = createCallbackFcn(app, @DEMcoregpathEditFieldValueChanged, true);
            app.DEMcoregpathEditField.HorizontalAlignment = 'right';
            app.DEMcoregpathEditField.FontName = 'Manrope';
            app.DEMcoregpathEditField.Position = [158 190 485 22];

            % Create DEMresamplingmethodDropDownLabel
            app.DEMresamplingmethodDropDownLabel = uilabel(app.SlavesProcessingTab);
            app.DEMresamplingmethodDropDownLabel.HorizontalAlignment = 'right';
            app.DEMresamplingmethodDropDownLabel.FontWeight = 'bold';
            app.DEMresamplingmethodDropDownLabel.Position = [41 150 146 22];
            app.DEMresamplingmethodDropDownLabel.Text = 'DEM resampling method';

            % Create DEMresamplingmethodDropDown
            app.DEMresamplingmethodDropDown = uidropdown(app.SlavesProcessingTab);
            app.DEMresamplingmethodDropDown.Items = {'NEAREST_NEIGHBOUR', 'BILINEAR_INTERPOLATION', 'CUBIC_CONVOLUTION', 'BISINC_5_POINT_INTERPOLATION', 'BISINC_11_POINT_INTERPOLATION', 'BISINC_21_POINT_INTERPOLATION', 'BICUBIC_INTERPOLATION'};
            app.DEMresamplingmethodDropDown.ValueChangedFcn = createCallbackFcn(app, @DEMresamplingmethodDropDownValueChanged, true);
            app.DEMresamplingmethodDropDown.FontName = 'Manrope';
            app.DEMresamplingmethodDropDown.FontWeight = 'bold';
            app.DEMresamplingmethodDropDown.BackgroundColor = [1 1 1];
            app.DEMresamplingmethodDropDown.Position = [197 150 201 22];
            app.DEMresamplingmethodDropDown.Value = 'NEAREST_NEIGHBOUR';

            % Create onlyforExternalDEMLabel
            app.onlyforExternalDEMLabel = uilabel(app.SlavesProcessingTab);
            app.onlyforExternalDEMLabel.FontName = 'Manrope';
            app.onlyforExternalDEMLabel.Position = [779 270 131 22];
            app.onlyforExternalDEMLabel.Text = '(only for External DEM)';

            % Create onlyforExternalDEMLabel_2
            app.onlyforExternalDEMLabel_2 = uilabel(app.SlavesProcessingTab);
            app.onlyforExternalDEMLabel_2.FontName = 'Manrope';
            app.onlyforExternalDEMLabel_2.Position = [779 190 131 22];
            app.onlyforExternalDEMLabel_2.Text = '(only for External DEM)';

            % Create FirststepDropDownLabel
            app.FirststepDropDownLabel = uilabel(app.SlavesProcessingTab);
            app.FirststepDropDownLabel.HorizontalAlignment = 'right';
            app.FirststepDropDownLabel.FontWeight = 'bold';
            app.FirststepDropDownLabel.Position = [41 396 59 22];
            app.FirststepDropDownLabel.Text = 'First step';

            % Create FirststepDropDown
            app.FirststepDropDown = uidropdown(app.SlavesProcessingTab);
            app.FirststepDropDown.Items = {'1', '2', '3', '4', '5', '6'};
            app.FirststepDropDown.ValueChangedFcn = createCallbackFcn(app, @FirststepDropDownValueChanged, true);
            app.FirststepDropDown.FontName = 'Manrope';
            app.FirststepDropDown.FontWeight = 'bold';
            app.FirststepDropDown.BackgroundColor = [1 1 1];
            app.FirststepDropDown.Position = [115 396 62 22];
            app.FirststepDropDown.Value = '1';

            % Create BrowseButton
            app.BrowseButton = uibutton(app.SlavesProcessingTab, 'push');
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);
            app.BrowseButton.BackgroundColor = [1 1 1];
            app.BrowseButton.FontName = 'Manrope';
            app.BrowseButton.Position = [662 269 100 23];
            app.BrowseButton.Text = 'Browse';

            % Create BrowseButton_2
            app.BrowseButton_2 = uibutton(app.SlavesProcessingTab, 'push');
            app.BrowseButton_2.ButtonPushedFcn = createCallbackFcn(app, @BrowseButton_2Pushed, true);
            app.BrowseButton_2.BackgroundColor = [1 1 1];
            app.BrowseButton_2.FontName = 'Manrope';
            app.BrowseButton_2.Position = [662 189 100 23];
            app.BrowseButton_2.Text = 'Browse';

            % Create CoherenceandLIATab
            app.CoherenceandLIATab = uitab(app.TabGroup);
            app.CoherenceandLIATab.Title = 'Coherence and LIA';
            app.CoherenceandLIATab.BackgroundColor = [1 1 1];

            % Create TerraincorrectedCoherenceandLIACheckBox
            app.TerraincorrectedCoherenceandLIACheckBox = uicheckbox(app.CoherenceandLIATab);
            app.TerraincorrectedCoherenceandLIACheckBox.ValueChangedFcn = createCallbackFcn(app, @TerraincorrectedCoherenceandLIACheckBoxValueChanged, true);
            app.TerraincorrectedCoherenceandLIACheckBox.Text = 'Terrain-corrected Coherence and LIA';
            app.TerraincorrectedCoherenceandLIACheckBox.FontName = 'Manrope';
            app.TerraincorrectedCoherenceandLIACheckBox.FontWeight = 'bold';
            app.TerraincorrectedCoherenceandLIACheckBox.Position = [48 437 238 22];
            app.TerraincorrectedCoherenceandLIACheckBox.Value = true;

            % Create EPSGcodeEditFieldLabel
            app.EPSGcodeEditFieldLabel = uilabel(app.CoherenceandLIATab);
            app.EPSGcodeEditFieldLabel.HorizontalAlignment = 'right';
            app.EPSGcodeEditFieldLabel.Position = [52 364 68 22];
            app.EPSGcodeEditFieldLabel.Text = 'EPSG code';

            % Create EPSGcodeEditField
            app.EPSGcodeEditField = uieditfield(app.CoherenceandLIATab, 'numeric');
            app.EPSGcodeEditField.Limits = [0 99999];
            app.EPSGcodeEditField.ValueDisplayFormat = '%.0f';
            app.EPSGcodeEditField.ValueChangedFcn = createCallbackFcn(app, @EPSGcodeEditFieldValueChanged, true);
            app.EPSGcodeEditField.FontName = 'Manrope';
            app.EPSGcodeEditField.Position = [135 364 100 22];
            app.EPSGcodeEditField.Value = 32633;

            % Create Label_3
            app.Label_3 = uilabel(app.CoherenceandLIATab);
            app.Label_3.FontName = 'Manrope';
            app.Label_3.FontWeight = 'bold';
            app.Label_3.Position = [48 395 670 22];
            app.Label_3.Text = 'EPSG code of the cartographic projection of the terrain-corrected .TIFF files (choose accordingly to the location)';

            % Create ComputationalResourcesTab
            app.ComputationalResourcesTab = uitab(app.TabGroup);
            app.ComputationalResourcesTab.Title = 'Computational Resources';
            app.ComputationalResourcesTab.BackgroundColor = [1 1 1];

            % Create FullpathtoSNAPgptfolderLabel
            app.FullpathtoSNAPgptfolderLabel = uilabel(app.ComputationalResourcesTab);
            app.FullpathtoSNAPgptfolderLabel.FontName = 'Manrope';
            app.FullpathtoSNAPgptfolderLabel.FontWeight = 'bold';
            app.FullpathtoSNAPgptfolderLabel.Position = [48 437 165 22];
            app.FullpathtoSNAPgptfolderLabel.Text = 'Full path to SNAP gpt folder';

            % Create PathEditFieldLabel
            app.PathEditFieldLabel = uilabel(app.ComputationalResourcesTab);
            app.PathEditFieldLabel.HorizontalAlignment = 'right';
            app.PathEditFieldLabel.Position = [54 406 30 22];
            app.PathEditFieldLabel.Text = 'Path';

            % Create PathEditField
            app.PathEditField = uieditfield(app.ComputationalResourcesTab, 'text');
            app.PathEditField.ValueChangedFcn = createCallbackFcn(app, @PathEditFieldValueChanged, true);
            app.PathEditField.FontName = 'Manrope';
            app.PathEditField.Position = [99 406 284 22];
            app.PathEditField.Value = 'C:\Program Files\snap\bin\gpt';

            % Create NumberofcorestobeusedintheprocessingLabel
            app.NumberofcorestobeusedintheprocessingLabel = uilabel(app.ComputationalResourcesTab);
            app.NumberofcorestobeusedintheprocessingLabel.FontName = 'Manrope';
            app.NumberofcorestobeusedintheprocessingLabel.FontWeight = 'bold';
            app.NumberofcorestobeusedintheprocessingLabel.Position = [45 359 269 22];
            app.NumberofcorestobeusedintheprocessingLabel.Text = 'Number of cores to be used in the processing';

            % Create CPUEditFieldLabel
            app.CPUEditFieldLabel = uilabel(app.ComputationalResourcesTab);
            app.CPUEditFieldLabel.HorizontalAlignment = 'right';
            app.CPUEditFieldLabel.Position = [57 325 30 22];
            app.CPUEditFieldLabel.Text = 'CPU';

            % Create CPUEditField
            app.CPUEditField = uieditfield(app.ComputationalResourcesTab, 'numeric');
            app.CPUEditField.Limits = [0 1000];
            app.CPUEditField.RoundFractionalValues = 'on';
            app.CPUEditField.ValueDisplayFormat = '%.0f';
            app.CPUEditField.ValueChangedFcn = createCallbackFcn(app, @CPUEditFieldValueChanged, true);
            app.CPUEditField.HorizontalAlignment = 'left';
            app.CPUEditField.FontName = 'Manrope';
            app.CPUEditField.Position = [102 325 73 22];
            app.CPUEditField.Value = 8;

            % Create RAMtobeusedintheprocessingGBformatnnGLabel
            app.RAMtobeusedintheprocessingGBformatnnGLabel = uilabel(app.ComputationalResourcesTab);
            app.RAMtobeusedintheprocessingGBformatnnGLabel.FontName = 'Manrope';
            app.RAMtobeusedintheprocessingGBformatnnGLabel.FontWeight = 'bold';
            app.RAMtobeusedintheprocessingGBformatnnGLabel.Position = [48 273 320 22];
            app.RAMtobeusedintheprocessingGBformatnnGLabel.Text = 'RAM to be used in the processing [GB] (format "nnG")';

            % Create CacheEditFieldLabel
            app.CacheEditFieldLabel = uilabel(app.ComputationalResourcesTab);
            app.CacheEditFieldLabel.HorizontalAlignment = 'right';
            app.CacheEditFieldLabel.Position = [55 242 40 22];
            app.CacheEditFieldLabel.Text = 'Cache';

            % Create CacheEditField
            app.CacheEditField = uieditfield(app.ComputationalResourcesTab, 'text');
            app.CacheEditField.ValueChangedFcn = createCallbackFcn(app, @CacheEditFieldValueChanged, true);
            app.CacheEditField.FontName = 'Manrope';
            app.CacheEditField.Position = [110 242 100 22];
            app.CacheEditField.Value = '26G';

            % Create ImagesTab_SEN
            app.ImagesTab_SEN = uitab(app.TabGroup);
            app.ImagesTab_SEN.Title = 'Images';
            app.ImagesTab_SEN.BackgroundColor = [1 1 1];

            % Create ImportedSENImagesTable
            app.ImportedSENImagesTable = uitable(app.ImagesTab_SEN);
            app.ImportedSENImagesTable.ColumnName = {'File'; 'Date'; 'Status'};
            app.ImportedSENImagesTable.RowName = {};
            app.ImportedSENImagesTable.Position = [262 66 675 325];

            % Create Label_SENImages
            app.Label_SENImages = uilabel(app.ImagesTab_SEN);
            app.Label_SENImages.HorizontalAlignment = 'center';
            app.Label_SENImages.WordWrap = 'on';
            app.Label_SENImages.FontName = 'Manrope';
            app.Label_SENImages.Position = [262 394 676 41];
            app.Label_SENImages.Text = 'Selected Sentinel-1 ZIP images will be copied/moved into PHASE_Preprocessing\slaves';

            % Create ImportSENButton
            app.ImportSENButton = uibutton(app.ImagesTab_SEN, 'push');
            app.ImportSENButton.ButtonPushedFcn = createCallbackFcn(app, @ImportSENButtonPushed, true);
            app.ImportSENButton.FontSize = 14;
            app.ImportSENButton.Position = [397 432 399 40];
            app.ImportSENButton.Text = 'Import Sentinel-1 images (.zip)';

            % Create OpenSENSlavesFolderButton
            app.OpenSENSlavesFolderButton = uibutton(app.ImagesTab_SEN, 'push');
            app.OpenSENSlavesFolderButton.ButtonPushedFcn = createCallbackFcn(app, @OpenSENSlavesFolderButtonPushed, true);
            app.OpenSENSlavesFolderButton.FontName = 'Manrope';
            app.OpenSENSlavesFolderButton.Position = [523 23 158 23];
            app.OpenSENSlavesFolderButton.Text = 'Open slaves folder';

            % Create SaveLoadTab
            app.SaveLoadTab = uitab(app.TabGroup);
            app.SaveLoadTab.Title = 'Save/Load';
            app.SaveLoadTab.BackgroundColor = [1 1 1];

            % Create SaveButton
            app.SaveButton = uibutton(app.SaveLoadTab, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.BackgroundColor = [1 1 1];
            app.SaveButton.FontName = 'Manrope';
            app.SaveButton.Position = [44 404 100 23];
            app.SaveButton.Text = 'Save';

            % Create SavetheconfiguredparametersfortheInSARpreprocessingLabel
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel = uilabel(app.SaveLoadTab);
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel.FontName = 'Manrope';
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel.FontWeight = 'bold';
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel.Position = [44 439 367 22];
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel.Text = 'Save the configured parameters for the InSAR pre-processing';

            % Create StatusLampLabel
            app.StatusLampLabel = uilabel(app.SaveLoadTab);
            app.StatusLampLabel.HorizontalAlignment = 'right';
            app.StatusLampLabel.Position = [171 404 39 22];
            app.StatusLampLabel.Text = 'Status';

            % Create StatusLamp
            app.StatusLamp = uilamp(app.SaveLoadTab);
            app.StatusLamp.Position = [225 404 20 20];
            app.StatusLamp.Color = [1 0 0];

            % Create SavetheconfiguredparametersfortheInSARpreprocessingLabel_3
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_3 = uilabel(app.SaveLoadTab);
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_3.FontName = 'Manrope';
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_3.FontWeight = 'bold';
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_3.Position = [45 328 502 22];
            app.SavetheconfiguredparametersfortheInSARpreprocessingLabel_3.Text = 'Load the parameters for the InSAR pre-processing from the already existing .mat file';

            % Create LoadButton
            app.LoadButton = uibutton(app.SaveLoadTab, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.BackgroundColor = [1 1 1];
            app.LoadButton.FontName = 'Manrope';
            app.LoadButton.Position = [43 293 100 23];
            app.LoadButton.Text = 'Load';

            % Create StatusLamp_3Label
            app.StatusLamp_3Label = uilabel(app.SaveLoadTab);
            app.StatusLamp_3Label.HorizontalAlignment = 'right';
            app.StatusLamp_3Label.Position = [171 294 39 22];
            app.StatusLamp_3Label.Text = 'Status';

            % Create StatusLamp_3
            app.StatusLamp_3 = uilamp(app.SaveLoadTab);
            app.StatusLamp_3.Position = [225 294 20 20];
            app.StatusLamp_3.Color = [1 0 0];

            % Create RunTab
            app.RunTab = uitab(app.TabGroup);
            app.RunTab.Title = 'Run';
            app.RunTab.BackgroundColor = [1 1 1];

            % Create Label_6
            app.Label_6 = uilabel(app.RunTab);
            app.Label_6.FontName = 'Manrope';
            app.Label_6.FontWeight = 'bold';
            app.Label_6.Position = [40 439 525 22];
            app.Label_6.Text = 'Press the button below to start the pre-processing of the SAR images for the PS analysis';

            % Create StartButton
            app.StartButton = uibutton(app.RunTab, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.FontName = 'Manrope';
            app.StartButton.Position = [40 406 100 23];
            app.StartButton.Text = 'Start';

            % Create Label_7
            app.Label_7 = uilabel(app.RunTab);
            app.Label_7.FontName = 'Manrope';
            app.Label_7.FontWeight = 'bold';
            app.Label_7.Position = [40 342 117 22];
            app.Label_7.Text = 'Execution outputs:';

            % Create PreprocessingstatusLampLabel
            app.PreprocessingstatusLampLabel = uilabel(app.RunTab);
            app.PreprocessingstatusLampLabel.HorizontalAlignment = 'right';
            app.PreprocessingstatusLampLabel.Position = [181 405 117 22];
            app.PreprocessingstatusLampLabel.Text = 'Preprocessing status';

            % Create PreprocessingstatusLamp
            app.PreprocessingstatusLamp = uilamp(app.RunTab);
            app.PreprocessingstatusLamp.Position = [313 405 20 20];
            app.PreprocessingstatusLamp.Color = [1 0 0];

            % Create MessagesTextArea_3Label
            app.MessagesTextArea_3Label = uilabel(app.RunTab);
            app.MessagesTextArea_3Label.HorizontalAlignment = 'right';
            app.MessagesTextArea_3Label.Position = [36 309 60 22];
            app.MessagesTextArea_3Label.Text = 'Messages';

            % Create MessagesTextArea
            app.MessagesTextArea = uitextarea(app.RunTab);
            app.MessagesTextArea.FontName = 'Manrope';
            app.MessagesTextArea.Position = [111 93 495 240];

            % Create StopButton
            app.StopButton = uibutton(app.RunTab, 'push');
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);
            app.StopButton.BackgroundColor = [1 0.8588 0.8588];
            app.StopButton.FontName = 'Manrope';
            app.StopButton.Position = [632 311 100 23];
            app.StopButton.Text = 'Stop';

            % Create ThecodewillstopattheendofthecurrentstepLabel_2
            app.ThecodewillstopattheendofthecurrentstepLabel_2 = uilabel(app.RunTab);
            app.ThecodewillstopattheendofthecurrentstepLabel_2.FontName = 'Manrope';
            app.ThecodewillstopattheendofthecurrentstepLabel_2.Position = [750 312 271 22];
            app.ThecodewillstopattheendofthecurrentstepLabel_2.Text = 'The code will stop at the end of the current step';

            % Create ConstellationSwitch
            app.ConstellationSwitch = uiswitch(app.UIFigure, 'slider');
            app.ConstellationSwitch.Items = {'Sentinel1', 'CosmoSkyMed'};
            app.ConstellationSwitch.ValueChangedFcn = createCallbackFcn(app, @ConstellationSwitchValueChanged, true);
            app.ConstellationSwitch.FontName = 'Manrope';
            app.ConstellationSwitch.FontSize = 14;
            app.ConstellationSwitch.Position = [562 598 45 20];
            app.ConstellationSwitch.Value = 'Sentinel1';

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.HorizontalAlignment = 'center';
            app.Label.FontName = 'Manrope';
            app.Label.FontSize = 14;
            app.Label.Position = [332 633 538 22];
            app.Label.Text = 'Move the switch to match the satellite constellation you want to process the data';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Position = [1032 557 154 105];
            app.Image.ImageSource = fullfile(pathToMLAPP, 'PHASE_logo.png');

            % Create Image2
            app.Image2 = uiimage(app.UIFigure);
            app.Image2.Position = [15 557 154 105];
            app.Image2.ImageSource = fullfile(pathToMLAPP, 'PHASE_mod1a.png');

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create Menu
            app.Menu = uimenu(app.ContextMenu);
            app.Menu.Text = 'Menu';

            % Create Menu2
            app.Menu2 = uimenu(app.ContextMenu);
            app.Menu2.Text = 'Menu2';

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
