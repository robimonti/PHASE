classdef LegacyAppAdapter < handle
    %LEGACYAPPADAPTER Text-based bridge for the proven processing callback.
    %
    % The first beta deliberately runs the same Start logic as the stable
    % MLAPP. These lightweight properties replace App Designer controls while
    % logging, status and TS-picker actions are delegated to the new app.

    properties
        Owner
        WorkDir
        UIFigure
        MessagesTextArea
        PreprocessingstatusLamp
        headingEditField
        lambdaEditField
        MasterdateDatePicker
        DateDatePicker
        DaysEditField
        TimeEditField
        subtr_tropoEditField
        tropo_methodEditField
        TRAINatmosphericcorrectionCheckBox
        weed_time_winEditField
        unwrap_time_winEditField
        scn_time_winEditField
        StaMPSfirststepDropDown
        StaMPSlaststepDropDown
        FilenameEditField
        master_date
        year_0
        month_0
        day_0
        train_flag
        TabGroup
        TSPointsTab
    end

    methods
        function obj = LegacyAppAdapter(owner, cfg)
            obj.Owner = owner;
            obj.WorkDir = owner.WorkDir;
            obj.UIFigure = owner.UIFigure;
            obj.MessagesTextArea = struct('Value', {{}});
            obj.PreprocessingstatusLamp = struct('Color', [1 0 0]);

            % Enable='on' means the MAT file is authoritative. The beta writes
            % auto-detected values to the visible configuration and requires an
            % explicit Save before Start, avoiding hidden runtime overrides.
            obj.headingEditField = field(cfg.heading);
            obj.lambdaEditField = field(cfg.lambda);
            obj.MasterdateDatePicker = field(cfg.master_date);
            obj.DateDatePicker = field(sprintf('%04d-%02d-%02d', ...
                cfg.year_0, cfg.month_0, cfg.day_0));
            obj.DaysEditField = field(cfg.time_span);
            obj.TimeEditField = field(cfg.utc_time);
            obj.subtr_tropoEditField = field(cfg.subtr_tropo);
            obj.tropo_methodEditField = field(cfg.tropo_method);
            obj.TRAINatmosphericcorrectionCheckBox = field(cfg.train_flag == 0);
            obj.weed_time_winEditField = field(cfg.weed_time_win);
            obj.unwrap_time_winEditField = field(cfg.unwrap_time_win);
            obj.scn_time_winEditField = field(cfg.scn_time_win);
            obj.StaMPSfirststepDropDown = field(char(string(cfg.stamps_first_step)));
            obj.StaMPSlaststepDropDown = field(char(string(cfg.stamps_last_step)));
            obj.FilenameEditField = field(cfg.export_name);
            obj.master_date = cfg.master_date;
            obj.year_0 = cfg.year_0;
            obj.month_0 = cfg.month_0;
            obj.day_0 = cfg.day_0;
            obj.train_flag = cfg.train_flag;
            obj.TabGroup = struct('SelectedTab', []);
            obj.TSPointsTab = [];
        end

        function log(obj, message)
            obj.Owner.appendLog(message);
        end

        function openTsPicker(obj)
            obj.Owner.openTsPicker();
        end
    end
end

function out = field(value)
out = struct('Value', value, 'Enable', 'on');
end
