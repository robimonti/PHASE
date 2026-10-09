function cleanup = renderTsPicker(workDir, host, valueType, exportName)
%RENDERTSPICKER Native map integrated into the StaMPS TS Points section.
% The complete export already exists; this view only exports chosen points.

matPath = fullfile(workDir,['ps_plot_ts_' valueType '.mat']);
series = load(matPath,'ph_mm','lonlat','day');
if size(series.ph_mm,1) ~= size(series.lonlat,1) || ...
        size(series.ph_mm,2) ~= numel(series.day) || ...
        any(~isfinite(series.day(:))) || numel(unique(series.day)) < 2
    error('PHASE_StaMPS_beta:pickerSeriesInvalid', ...
        'The StaMPS time-series MAT file has inconsistent points or dates.');
end
if ~license('test','Map_Toolbox')
    error('PHASE_StaMPS_beta:pickerMapToolbox', ...
        'Mapping Toolbox is required for the interactive TS Points map.');
end

years = (double(series.day(:))-double(series.day(1)))/365.25;
rate = nan(size(series.lonlat,1),1);
centeredYears = years-mean(years);
completeRows = all(isfinite(series.ph_mm),2);
rate(completeRows) = double(series.ph_mm(completeRows,:))*centeredYears / ...
    sum(centeredYears.^2);
for k = find(~completeRows)'
    values = double(series.ph_mm(k,:))';
    valid = isfinite(values) & isfinite(years);
    if nnz(valid) >= 2
        fit = polyfit(years(valid),values(valid),1);
        rate(k) = fit(1);
    end
end

delete(host.Children);
blue = [0.208 0.396 0.812];
navy = [0.078 0.149 0.275];
muted = [0.39 0.47 0.61];
surface = [0.965 0.978 0.997];
host.BackgroundColor = surface;
previousScrollable = host.Scrollable;
host.Scrollable = 'on';
topfig = ancestor(host,'figure');
previousClick = topfig.WindowButtonDownFcn;
previousResize = host.SizeChangedFcn;
previousAutoResize = host.AutoResizeChildren;
pickMode = false;
selectedRow = [];
markers = gobjects(0);

% A grid that only fills the viewport never overflows, so Scrollable alone
% has no effect. Give the picker a minimum canvas height and resize it with
% the viewport; short remote-desktop windows can then scroll vertically.
canvas = uipanel(host,'BorderType','none','Units','pixels', ...
    'BackgroundColor',surface);
host.AutoResizeChildren = 'off';
host.SizeChangedFcn = @resizeCanvas;
resizeCanvas([],[]);
root = uigridlayout(canvas,[2 2]);
root.RowHeight = {77,'1x'};
root.ColumnWidth = {'1x',350};
root.Padding = [23 17 23 22];
root.RowSpacing = 16;
root.ColumnSpacing = 16;
root.BackgroundColor = surface;

header = uigridlayout(root,[2 1]);
header.Layout.Row = 1;
header.Layout.Column = [1 2];
header.RowHeight = {31,23};
header.Padding = [0 0 0 0];
header.RowSpacing = 1;
header.BackgroundColor = surface;
uilabel(header,'Text','TS Points','FontSize',23,'FontWeight','bold', ...
    'FontColor',navy);
uilabel(header,'Text', ...
    'Optional point export · The complete time series is already in project results.', ...
    'FontSize',11,'FontColor',muted);

mapCard = uipanel(root,'BorderType','line', ...
    'BackgroundColor',[1 1 1],'HighlightColor',[0.83 0.88 0.96]);
mapCard.Layout.Row = 2;
mapCard.Layout.Column = 1;
mapGrid = uigridlayout(mapCard,[2 1]);
mapGrid.RowHeight = {39,'1x'};
mapGrid.Padding = [13 12 13 13];
mapGrid.RowSpacing = 5;
mapGrid.BackgroundColor = [1 1 1];
mapHead = uigridlayout(mapGrid,[1 2]);
mapHead.Layout.Row = 1;
mapHead.ColumnWidth = {'1x','fit'};
mapHead.Padding = [4 0 2 0];
mapHead.BackgroundColor = [1 1 1];
uilabel(mapHead,'Text','PERSISTENT SCATTERERS', ...
    'FontSize',11,'FontWeight','bold','FontColor',navy);
mapHint = uilabel(mapHead,'Text','Navigate: pan and zoom', ...
    'FontSize',10,'FontColor',muted);
mapHint.Layout.Column = 2;

% A panel between geoaxes and the grid preserves map pan/zoom in R2026a.
mapHost = uipanel(mapGrid,'BorderType','none','BackgroundColor',[1 1 1]);
mapHost.Layout.Row = 2;
ax = geoaxes(mapHost);
geoscatter(ax,series.lonlat(:,2),series.lonlat(:,1),18,rate,'filled');
try
    geobasemap(ax,'satellite');
catch
    geobasemap(ax,'streets-light');
end
colourScale = colorbar(ax);
colourScale.Label.String = 'LOS rate [mm/yr]';
colormap(ax,turbo);
hold(ax,'on');
defaultInteractions = ax.Interactions;

side = uipanel(root,'BorderType','line', ...
    'BackgroundColor',[1 1 1],'HighlightColor',[0.83 0.88 0.96]);
side.Layout.Row = 2;
side.Layout.Column = 2;
controls = uigridlayout(side,[9 1]);
controls.RowHeight = {42,'1x',44,42,59,46,46,51,55};
controls.Padding = [15 13 15 13];
controls.RowSpacing = 8;
controls.BackgroundColor = [1 1 1];

uilabel(controls,'Text','SELECTED POINTS', ...
    'FontSize',11,'FontWeight','bold','FontColor',navy);
points = uitable(controls,'Data',cell(0,4), ...
    'ColumnName',{'ID','Longitude','Latitude','Radius m'}, ...
    'ColumnFormat',{'char','numeric','numeric','numeric'}, ...
    'ColumnEditable',[true true true true], ...
    'ColumnWidth',{64,75,75,70}, ...
    'CellSelectionCallback',@selectRow);
points.Layout.Row = 2;

pick = actionButton(controls,'◎  Select a scatterer on the map', ...
    [0.92 0.95 1],blue,@togglePick);
pick.Layout.Row = 3;
free = actionButton(controls,'＋  Place a point anywhere on the map', ...
    [0.96 0.97 0.99],navy,@addFreePoint);
free.Layout.Row = 4;

coordinates = uigridlayout(controls,[2 3]);
coordinates.Layout.Row = 5;
coordinates.RowHeight = {16,32};
coordinates.ColumnWidth = {'1x','1x','1x'};
coordinates.Padding = [0 0 0 0];
coordinates.ColumnSpacing = 5;
for labelIndex = 1:3
    labels = {'LONGITUDE','LATITUDE','RADIUS M'};
    label = uilabel(coordinates,'Text',labels{labelIndex}, ...
        'FontSize',9,'FontWeight','bold','FontColor',muted);
    label.Layout.Row = 1;
    label.Layout.Column = labelIndex;
end
lonInput = uieditfield(coordinates,'numeric','Value',0, ...
    'Limits',[-180 180]); lonInput.Layout.Row = 2; lonInput.Layout.Column = 1;
latInput = uieditfield(coordinates,'numeric','Value',0, ...
    'Limits',[-90 90]); latInput.Layout.Row = 2; latInput.Layout.Column = 2;
radiusInput = uieditfield(coordinates,'numeric','Value',100, ...
    'Limits',[1 1e6]); radiusInput.Layout.Row = 2; radiusInput.Layout.Column = 3;

manual = actionButton(controls,'Add typed coordinates', ...
    [0.92 0.95 1],blue,@addManual);
manual.Layout.Row = 6;
remove = actionButton(controls,'Remove highlighted point', ...
    [0.96 0.97 0.99],navy,@removeSelected);
remove.Layout.Row = 7;

fileActions = uigridlayout(controls,[1 2]);
fileActions.Layout.Row = 8;
fileActions.ColumnWidth = {'1x','1x'};
fileActions.Padding = [0 0 0 0];
loadButton = actionButton(fileActions,'Import points CSV', ...
    [0.92 0.95 1],blue,@loadList);
loadButton.Layout.Column = 1;
saveButton = actionButton(fileActions,'Save points CSV', ...
    [0.92 0.95 1],blue,@saveList);
saveButton.Layout.Column = 2;

bottom = uigridlayout(controls,[2 1]);
bottom.Layout.Row = 9;
bottom.RowHeight = {32,18};
bottom.Padding = [0 0 0 0];
bottom.RowSpacing = 3;
exportButton = actionButton(bottom,'Export selected points', ...
    blue,[1 1 1],@exportSelected);
exportButton.Layout.Row = 1;
status = uilabel(bottom,'Text','Ready · select a point or load a list.', ...
    'FontSize',9,'FontColor',muted);
status.Layout.Row = 2;

autoload = fullfile(workDir,'aoi_points.csv');
if isfile(autoload)
    try
        readList(autoload);
        status.Text = 'Previous point list loaded.';
    catch ME
        status.Text = ['Point list not loaded: ' ME.message];
    end
end
cleanup = @restoreMap;

    function button = actionButton(parent,textValue,background,foreground,callback)
        htmlSource = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
            'phase_stamps_beta_ui','rounded_button.html');
        button = uihtml(parent,'HTMLSource',htmlSource, ...
            'Data',struct('label',textValue, ...
                'background',cssColor(background), ...
                'foreground',cssColor(foreground), ...
                'enabled',true,'clicked',0), ...
            'DataChangedFcn',callback);
    end

    function value = cssColor(rgb)
        channels = round(255*rgb);
        value = sprintf('#%02X%02X%02X',channels(1),channels(2),channels(3));
    end

    function setButtonState(button,field,value)
        data = button.Data;
        data.(field) = value;
        button.Data = data;
    end

    function restoreMap()
        try
            topfig.WindowButtonDownFcn = previousClick;
            host.SizeChangedFcn = previousResize;
            host.AutoResizeChildren = previousAutoResize;
            host.Scrollable = previousScrollable;
            if isvalid(ax), ax.Interactions = defaultInteractions; end
        catch
        end
    end

    function resizeCanvas(~,~)
        if ~isvalid(host) || ~isvalid(canvas), return; end
        viewport = host.Position;
        canvas.Position = [0 0 max(760,viewport(3)-14) ...
            max(710,viewport(4)-14)];
    end

    function togglePick(~,~)
        pickMode = ~pickMode;
        if pickMode
            ax.Interactions = zoomInteraction;
            topfig.WindowButtonDownFcn = @onMapClick;
            setButtonState(pick,'label','●  Picking scatterers · click a dot');
            setButtonState(pick,'background',cssColor([0.83 0.91 1]));
            mapHint.Text = 'Pick mode · click a PS';
        else
            topfig.WindowButtonDownFcn = previousClick;
            ax.Interactions = defaultInteractions;
            setButtonState(pick,'label','◎  Select a scatterer on the map');
            setButtonState(pick,'background',cssColor([0.92 0.95 1]));
            mapHint.Text = 'Navigate: pan and zoom';
        end
    end

    function onMapClick(~,~)
        if ~pickMode, return; end
        cursor = topfig.CurrentPoint;
        box = getpixelposition(ax,true);
        if cursor(1) < box(1) || cursor(1) > box(1)+box(3) || ...
                cursor(2) < box(2) || cursor(2) > box(2)+box(4)
            return
        end
        point = ax.CurrentPoint; % geoaxes order is latitude, longitude.
        lon = point(1,2);
        lat = point(1,1);
        dx = (double(series.lonlat(:,1))-lon).*111320.*cosd(lat);
        dy = (double(series.lonlat(:,2))-lat).*111320;
        [~,nearest] = min(dx.^2+dy.^2);
        addPoint(nextId(),series.lonlat(nearest,1), ...
            series.lonlat(nearest,2),radiusInput.Value);
        status.Text = sprintf('Added PS #%d.',nearest);
    end

    function addFreePoint(~,~)
        if pickMode, togglePick([],[]); end
        status.Text = 'Click a position on the map.';
        try
            roi = drawpoint(ax,'Color',blue);
            if isempty(roi.Position), delete(roi); return; end
            location = roi.Position;
            delete(roi);
            addPoint(nextId(),location(2),location(1),radiusInput.Value);
            status.Text = 'Free point added.';
        catch ME
            status.Text = ['Map selection cancelled: ' ME.message];
        end
    end

    function addManual(~,~)
        addPoint(nextId(),lonInput.Value,latInput.Value,radiusInput.Value);
        status.Text = 'Coordinates added to the selection.';
    end

    function addPoint(id,lon,lat,radius)
        row = {char(string(id)),double(lon),double(lat),double(radius)};
        points.Data = [points.Data; row];
        [latLimits,lonLimits] = geolimits(ax);
        marker = geoplot(ax,lat,lon,'p', ...
            'MarkerSize',12,'LineWidth',1.8,'Color',[0.89 0.17 0.40]);
        marker.HitTest = 'off';
        marker.PickableParts = 'none';
        markers(end+1) = marker;
        geolimits(ax,latLimits,lonLimits);
    end

    function id = nextId()
        ids = string(points.Data(:,1));
        n = 1;
        while any(ids == string(sprintf('P%02d',n))), n = n+1; end
        id = sprintf('P%02d',n);
    end

    function selectRow(~,event)
        selectedRow = [];
        if ~isempty(event.Indices), selectedRow = event.Indices(1,1); end
    end

    function removeSelected(~,~)
        if isempty(selectedRow) || selectedRow > size(points.Data,1)
            status.Text = 'Select a table row first.';
            return
        end
        points.Data(selectedRow,:) = [];
        if selectedRow <= numel(markers)
            if isgraphics(markers(selectedRow)), delete(markers(selectedRow)); end
            markers(selectedRow) = [];
        end
        selectedRow = [];
        status.Text = 'Point removed.';
    end

    function loadList(~,~)
        [file,folder] = uigetfile({'*.csv','CSV files'}, ...
            'Load point list',workDir);
        if isequal(file,0), return; end
        try
            readList(fullfile(folder,file));
            status.Text = sprintf('Loaded %d points.',size(points.Data,1));
        catch ME
            uialert(topfig,ME.message,'Cannot load point list');
        end
    end

    function readList(pathValue)
        tableValue = readtable(pathValue,'TextType','string');
        required = {'id','lon','lat'};
        if ~all(ismember(required,tableValue.Properties.VariableNames))
            error('PHASE_StaMPS_beta:pickerCsv', ...
                'The CSV must contain id, lon and lat columns.');
        end
        if ~ismember('radius_m',tableValue.Properties.VariableNames)
            tableValue.radius_m = repmat(radiusInput.Value,height(tableValue),1);
        end
        if any(~isfinite(tableValue.lon) | abs(tableValue.lon)>180 | ...
               ~isfinite(tableValue.lat) | abs(tableValue.lat)>90 | ...
               ~isfinite(tableValue.radius_m) | tableValue.radius_m<=0)
            error('PHASE_StaMPS_beta:pickerCsv', ...
                'The point list contains invalid coordinates or radii.');
        end
        delete(markers(isgraphics(markers)));
        markers = gobjects(0);
        points.Data = cell(0,4);
        for rowIndex = 1:height(tableValue)
            addPoint(tableValue.id(rowIndex),tableValue.lon(rowIndex), ...
                tableValue.lat(rowIndex),tableValue.radius_m(rowIndex));
        end
        selectedRow = [];
    end

    function saveList(~,~)
        if isempty(points.Data)
            status.Text = 'Select at least one point first.';
            return
        end
        [file,folder] = uiputfile({'*.csv','CSV files'}, ...
            'Save point list',fullfile(workDir,'aoi_points.csv'));
        if isequal(file,0), return; end
        writetable(selectionTable(),fullfile(folder,file));
        status.Text = 'Point list saved.';
    end

    function result = selectionTable()
        rows = points.Data;
        result = table(string(rows(:,1)),cell2mat(rows(:,2)), ...
            cell2mat(rows(:,3)),cell2mat(rows(:,4)), ...
            'VariableNames',{'id','lon','lat','radius_m'});
    end

    function exportSelected(~,~)
        if isempty(points.Data)
            status.Text = 'Select at least one point first.';
            return
        end
        setButtonState(exportButton,'enabled',false);
        status.Text = 'Exporting selected time series…';
        drawnow;
        temporaryCsv = [tempname '.csv'];
        temporaryCleanup = onCleanup(@() deleteIfPresent(temporaryCsv)); %#ok<NASGU>
        try
            writetable(selectionTable(),temporaryCsv);
            ts_export_batch(matPath,temporaryCsv, ...
                fullfile(workDir,'EXPORT'),radiusInput.Value,exportName);
            published = phase_stamps_beta.publishExports(workDir);
            if isempty(published)
                status.Text = sprintf('%d time series exported to EXPORT.',size(points.Data,1));
            else
                status.Text = sprintf('%d time series exported to project results.',size(points.Data,1));
            end
        catch ME
            status.Text = ['Export failed: ' ME.message];
        end
        setButtonState(exportButton,'enabled',true);
    end
end

function deleteIfPresent(pathValue)
if isfile(pathValue), delete(pathValue); end
end
