function report = exportGridDiagnostics(workDir, cfg, logFcn)
%EXPORTGRIDDIAGNOSTICS Export radar sampling grid and PS selection history.
%
% The grid is the SLC/interferometric sampling grid, not the physical SAR
% resolution. SNAP longitude/latitude rasters are preferred. Candidate
% interpolation is used only when those source rasters are unavailable.

if nargin < 3 || isempty(logFcn)
    logFcn = @(message) fprintf('%s\n', message);
end
workDir = char(java.io.File(workDir).getCanonicalPath());
exportDir = fullfile(workDir, 'EXPORT');
if ~isfolder(exportDir), mkdir(exportDir); end

[width, imageLength] = readDimensions(workDir);
threshold = readAmplitudeThreshold(workDir, cfg.amplitude_threshold);
patches = dir(fullfile(workDir, 'PATCH_*'));
patches = patches([patches.isdir]);
patches = patches(~ismember({patches.name}, {'.','..'}));
if isempty(patches)
    error('PHASE_StaMPS:gridDiagnosticNoPatches', ...
        'No PATCH_* folders were found in %s.', workDir);
end

records = repmat(emptyRecord(), 0, 1);
for k = 1:numel(patches)
    patchDir = fullfile(patches(k).folder, patches(k).name);
    record = loadPatchRecord(patchDir, patches(k).name, width, imageLength);
    records(end+1) = record; %#ok<AGROW>
end

[candidateTable, candidates] = concatenateRecords(records, threshold);
if isempty(candidates.lon)
    error('PHASE_StaMPS:gridDiagnosticNoCandidates', ...
        'The PATCH_* folders contain no initial PS candidates.');
end

[geo, geolocationMode, geoDetails] = loadGeolocationGrid( ...
    workDir, cfg, width, imageLength, candidates);
candidateTable.GeolocationMode(:) = string(geolocationMode);

runName = safeRunName(cfg.export_name, workDir);
candidateTable.RunName(:) = string(runName);
csvPath = fullfile(exportDir, ...
    sprintf('DiagnosticaGriglia_Stadi_%s.csv', runName));
writetable(candidateTable, csvPath);

maxGridLines = 150;
bufferMetres = 25;
[lineStepAz, lineStepRg] = gridSteps(imageLength, width, maxGridLines);
extent = bufferedExtent(candidates.lon, candidates.lat, bufferMetres);

daPath = fullfile(exportDir, sprintf('DiagnosticaGriglia_DA_%s.png', runName));
stagePath = fullfile(exportDir, sprintf('DiagnosticaGriglia_Stadi_%s.png', runName));
createDaFigure(daPath, geo, candidates, extent, lineStepAz, lineStepRg, ...
    threshold, runName, geolocationMode);
createStageFigure(stagePath, geo, candidates, extent, lineStepAz, lineStepRg, ...
    threshold, runName, geolocationMode);

metadata = struct();
metadata.runName = runName;
metadata.createdUtc = char(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
metadata.gridMeaning = ['SLC/interferometric sampling grid; it does not ' ...
    'represent the physical SAR resolution or point-spread function.'];
metadata.indexConvention = '[internal ID, zero-based azimuth pixel, zero-based range pixel]';
metadata.amplitudeDispersionThreshold = threshold;
metadata.thresholdSource = thresholdSource(workDir);
metadata.geolocationMode = geolocationMode;
metadata.geolocationDetails = geoDetails;
metadata.width = width;
metadata.length = imageLength;
metadata.patchCount = numel(records);
metadata.candidateCount = height(candidateTable);
metadata.rejectedByPsSelect = nnz(candidates.stage == 1);
metadata.rejectedByWeeding = nnz(candidates.stage == 2);
metadata.survivesPatchWeeding = nnz(candidates.stage == 3);
metadata.mergeResampleSizeMetres = cfg.merge_resample_size;
metadata.gridLineStepAzimuthPixels = lineStepAz;
metadata.gridLineStepRangePixels = lineStepRg;
metadata.maxGridLinesPerDirection = maxGridLines;
metadata.bufferMetres = bufferMetres;
metadata.outputs = struct('amplitudeDispersionFigure', daPath, ...
    'processingStageFigure', stagePath, 'candidateCsv', csvPath);
metadataPath = fullfile(exportDir, ...
    sprintf('DiagnosticaGriglia_Metadata_%s.json', runName));
writeJson(metadataPath, metadata);

logFcn(sprintf(['Grid diagnostic used %s geolocation; grid lines every ' ...
    '%d azimuth and %d range cells.'], geolocationMode, lineStepAz, lineStepRg));
logFcn(['Grid diagnostic figures: ' daPath ' ; ' stagePath]);
logFcn(['Grid diagnostic candidate table: ' csvPath]);

report = struct('ok', true, 'candidateCount', height(candidateTable), ...
    'patchCount', numel(records), 'geolocationMode', geolocationMode, ...
    'daFigure', daPath, 'stageFigure', stagePath, 'csv', csvPath, ...
    'metadata', metadataPath);
end

function record = emptyRecord()
record = struct('patch', '', 'internalId', [], 'azimuth', [], 'range', [], ...
    'lon', [], 'lat', [], 'da', [], 'passesSelect', [], 'survivesWeeding', []);
end

function record = loadPatchRecord(patchDir, patchName, width, imageLength)
required = {'ps1.mat','da1.mat','select1.mat','weed1.mat'};
for k = 1:numel(required)
    path = fullfile(patchDir, required{k});
    if exist(path, 'file') ~= 2
        error('PHASE_StaMPS:gridDiagnosticMissingStage', ...
            ['%s is missing. The processing-stage diagnostic requires ', ...
             'StaMPS steps 1 through 4 to have completed.'], path);
    end
end

ps = load(fullfile(patchDir, 'ps1.mat'), 'ij', 'lonlat', 'n_ps');
da = load(fullfile(patchDir, 'da1.mat'), 'D_A');
selection = load(fullfile(patchDir, 'select1.mat'), 'ix', 'keep_ix');
weed = load(fullfile(patchDir, 'weed1.mat'), 'ix_weed');

n = size(ps.ij, 1);
if isfield(ps, 'n_ps') && double(ps.n_ps) ~= n
    error('PHASE_StaMPS:gridDiagnosticPsCount', ...
        '%s/ps1.mat has inconsistent n_ps and ij dimensions.', patchName);
end
if size(ps.ij,2) < 3 || size(ps.lonlat,1) ~= n || size(ps.lonlat,2) < 2
    error('PHASE_StaMPS:gridDiagnosticPsShape', ...
        '%s/ps1.mat has invalid ij or lonlat dimensions.', patchName);
end
if numel(da.D_A) ~= n
    error('PHASE_StaMPS:gridDiagnosticDaCount', ...
        '%s/da1.mat contains %d D_A values for %d candidates.', ...
        patchName, numel(da.D_A), n);
end

selectionIx = double(selection.ix(:));
keep = logical(selection.keep_ix(:));
if numel(selectionIx) ~= numel(keep) || any(selectionIx < 1 | selectionIx > n)
    error('PHASE_StaMPS:gridDiagnosticSelectionIndex', ...
        '%s/select1.mat contains inconsistent candidate indices.', patchName);
end
selected = selectionIx(keep);
weedKeep = logical(weed.ix_weed(:));
if numel(weedKeep) ~= numel(selected)
    error('PHASE_StaMPS:gridDiagnosticWeedCount', ...
        ['%s/weed1.mat has %d ix_weed values, but %d candidates survived ', ...
         'ps_select.'], patchName, numel(weedKeep), numel(selected));
end

azimuth = double(ps.ij(:,2));
range = double(ps.ij(:,3));
if any(azimuth < 0 | azimuth >= imageLength | range < 0 | range >= width)
    error('PHASE_StaMPS:gridDiagnosticRadarIndex', ...
        ['%s contains radar indices outside the zero-based %d-by-%d ', ...
         'azimuth/range grid.'], patchName, imageLength, width);
end
passesSelect = false(n,1);
passesSelect(selected) = true;
survivesWeeding = false(n,1);
survivesWeeding(selected(weedKeep)) = true;

record = struct('patch', patchName, 'internalId', double(ps.ij(:,1)), ...
    'azimuth', azimuth, 'range', range, ...
    'lon', double(ps.lonlat(:,1)), 'lat', double(ps.lonlat(:,2)), ...
    'da', double(da.D_A(:)), 'passesSelect', passesSelect, ...
    'survivesWeeding', survivesWeeding);
end

function [out, values] = concatenateRecords(records, threshold)
patch = strings(0,1); patchIndex = zeros(0,1); internalId = zeros(0,1);
az = zeros(0,1); rg = zeros(0,1); lon = zeros(0,1); lat = zeros(0,1);
da = zeros(0,1); passes = false(0,1); survives = false(0,1);
for k = 1:numel(records)
    n = numel(records(k).da);
    patch = [patch; repmat(string(records(k).patch), n, 1)]; %#ok<AGROW>
    patchIndex = [patchIndex; (1:n)']; %#ok<AGROW>
    internalId = [internalId; records(k).internalId(:)]; %#ok<AGROW>
    az = [az; records(k).azimuth(:)]; %#ok<AGROW>
    rg = [rg; records(k).range(:)]; %#ok<AGROW>
    lon = [lon; records(k).lon(:)]; %#ok<AGROW>
    lat = [lat; records(k).lat(:)]; %#ok<AGROW>
    da = [da; records(k).da(:)]; %#ok<AGROW>
    passes = [passes; records(k).passesSelect(:)]; %#ok<AGROW>
    survives = [survives; records(k).survivesWeeding(:)]; %#ok<AGROW>
end
stage = ones(size(da));
stage(passes) = 2;
stage(survives) = 3;
labels = strings(size(stage));
labels(stage == 1) = "Eliminato in ps_select";
labels(stage == 2) = "Eliminato al weeding";
labels(stage == 3) = "Sopravvive al weeding (patch-local ps2)";
runNames = repmat("", size(stage));
thresholds = repmat(threshold, size(stage));
geoModes = repmat("pending", size(stage));
out = table(runNames, patch, patchIndex, internalId, az, rg, lon, lat, da, ...
    thresholds, passes, survives, stage, labels, geoModes, ...
    'VariableNames', {'RunName','PatchName','CandidateIndex','InternalId', ...
    'AzimuthPixel','RangePixel','Longitude','Latitude','D_A','Threshold_DA', ...
    'PassesPsSelect','SurvivesWeeding','LastProcessingStageCode', ...
    'LastProcessingStage','GeolocationMode'});
values = struct('patch',patch,'azimuth',az,'range',rg,'lon',lon,'lat',lat, ...
    'da',da,'passesSelect',passes,'survivesWeeding',survives,'stage',stage);
end

function [width, imageLength] = readDimensions(workDir)
width = readPositiveInteger(fullfile(workDir, 'width.txt'), 'width');
imageLength = readPositiveInteger(fullfile(workDir, 'len.txt'), 'length');
end

function value = readPositiveInteger(path, label)
if exist(path, 'file') ~= 2
    error('PHASE_StaMPS:gridDiagnosticDimensionMissing', ...
        'Radar-grid %s file is missing: %s', label, path);
end
value = str2double(strtrim(fileread(path)));
if ~isfinite(value) || value < 1 || value ~= round(value)
    error('PHASE_StaMPS:gridDiagnosticDimensionInvalid', ...
        'Radar-grid %s is invalid in %s.', label, path);
end
end

function threshold = readAmplitudeThreshold(workDir, fallback)
threshold = fallback;
path = fullfile(workDir, 'selpsc.in');
if exist(path, 'file') == 2
    fid = fopen(path, 'r');
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    line = fgetl(fid);
    parsed = str2double(strtrim(line));
    if isfinite(parsed) && parsed >= 0, threshold = parsed; end
end
end

function source = thresholdSource(workDir)
if exist(fullfile(workDir,'selpsc.in'),'file') == 2
    source = 'selpsc.in';
else
    source = 'input_StaMPS.mat amplitude_threshold';
end
end

function [geo, mode, details] = loadGeolocationGrid(workDir, cfg, width, imageLength, candidates)
[lonPath, latPath] = resolveGeolocationPaths(workDir, cfg);
if ~isempty(lonPath) && ~isempty(latPath)
    [lonGrid, latGrid, byteOrder, residual] = readAndValidateGeoRasters( ...
        lonPath, latPath, width, imageLength, candidates);
    geo = struct('mode','grid','lon',lonGrid,'lat',latGrid, ...
        'lonF',griddedInterpolant({1:imageLength,1:width},lonGrid,'linear','linear'), ...
        'latF',griddedInterpolant({1:imageLength,1:width},latGrid,'linear','linear'), ...
        'width',width,'length',imageLength);
    mode = 'original_grid';
    details = struct('longitudeRaster',lonPath,'latitudeRaster',latPath, ...
        'byteOrder',byteOrder,'candidateCentreMedianResidualDegrees',residual);
    return
end

[uniqueIj, ~, group] = unique([candidates.azimuth candidates.range], 'rows');
lonMean = accumarray(group, candidates.lon, [], @mean);
latMean = accumarray(group, candidates.lat, [], @mean);
if size(uniqueIj,1) < 3 || rank([uniqueIj ones(size(uniqueIj,1),1)]) < 3
    error('PHASE_StaMPS:gridDiagnosticGeolocationMissing', ...
        ['SNAP .lon/.lat rasters were not found and the saved candidates do ', ...
         'not span a two-dimensional area for the documented fallback.']);
end
geo = struct('mode','interpolant', ...
    'lon',scatteredInterpolant(uniqueIj(:,2),uniqueIj(:,1),lonMean,'natural','linear'), ...
    'lat',scatteredInterpolant(uniqueIj(:,2),uniqueIj(:,1),latMean,'natural','linear'), ...
    'width',width,'length',imageLength);
mode = 'interpolated_from_ps1';
details = struct('longitudeRaster','','latitudeRaster','', ...
    'byteOrder','','candidateCentreMedianResidualDegrees',0);
end

function [lonPath, latPath] = resolveGeolocationPaths(workDir, cfg)
lonPath = ''; latPath = '';
inputPath = fullfile(workDir, 'psclonlat.in');
if exist(inputPath, 'file') == 2
    lines = regexp(fileread(inputPath), '\r?\n', 'split');
    lines = lines(~cellfun(@isempty, strtrim(lines)));
    if numel(lines) >= 3
        candidateLon = stripQuotes(strtrim(lines{2}));
        candidateLat = stripQuotes(strtrim(lines{3}));
        if exist(candidateLon,'file') == 2 && exist(candidateLat,'file') == 2
            lonPath = candidateLon; latPath = candidateLat; return
        end
    end
end

master = char(string(cfg.master_date));
root = char(string(cfg.project_path));
projectPreprocessing = fullfile(root,'processing','preprocessing');
if ~isempty(phase_project.findRoot(root))
    projectPreprocessing = phase_project.paths(root).preprocessing;
end
folders = { ...
    fullfile(projectPreprocessing,['INSAR_' master],'geo'), ...
    fullfile(root,'PHASE_Preprocessing',['INSAR_' master],'geo'), ...
    fullfile(root,['INSAR_' master],'geo'), ...
    fullfile(root,'engine','PHASE_Preprocessing',['INSAR_' master],'geo'), ...
    fullfile(fileparts(workDir),'PHASE_Preprocessing',['INSAR_' master],'geo')};
for k = 1:numel(folders)
    candidateLon = fullfile(folders{k}, [master '.lon']);
    candidateLat = fullfile(folders{k}, [master '.lat']);
    if exist(candidateLon,'file') == 2 && exist(candidateLat,'file') == 2
        lonPath = candidateLon; latPath = candidateLat; return
    end
end
end

function value = stripQuotes(value)
value = regexprep(value, '^["'']|["'']$', '');
end

function [lonGrid, latGrid, byteOrder, residual] = readAndValidateGeoRasters( ...
        lonPath, latPath, width, imageLength, candidates)
expectedBytes = width * imageLength * 4;
lonInfo = dir(lonPath); latInfo = dir(latPath);
if lonInfo.bytes ~= expectedBytes || latInfo.bytes ~= expectedBytes
    error('PHASE_StaMPS:gridDiagnosticGeoSize', ...
        ['SNAP geolocation rasters must contain exactly width*length float32 ', ...
         'values (%d bytes expected; lon=%d, lat=%d).'], ...
        expectedBytes, lonInfo.bytes, latInfo.bytes);
end

orders = {'ieee-le','ieee-be'};
bestScore = Inf; bestLon = []; bestLat = []; byteOrder = '';
sample = unique(round(linspace(1,numel(candidates.lon),min(500,numel(candidates.lon)))));
rows = candidates.azimuth(sample) + 1;
cols = candidates.range(sample) + 1;
validReference = isfinite(candidates.lon(sample)) & isfinite(candidates.lat(sample));
for k = 1:numel(orders)
    lon = readFloatRaster(lonPath,width,imageLength,orders{k});
    lat = readFloatRaster(latPath,width,imageLength,orders{k});
    indices = sub2ind([imageLength width],rows,cols);
    delta = hypot(lon(indices)-candidates.lon(sample), ...
        lat(indices)-candidates.lat(sample));
    delta = delta(validReference & isfinite(delta));
    if isempty(delta), score = Inf; else, score = median(delta); end
    if score < bestScore
        bestScore = score; bestLon = lon; bestLat = lat; byteOrder = orders{k};
    end
end
if isempty(bestLon) || bestScore > 0.01
    error('PHASE_StaMPS:gridDiagnosticGeoAlignment', ...
        ['SNAP .lon/.lat rasters could not be aligned with ps1 candidates. ', ...
         'Best median centre residual was %.6g degrees.'], bestScore);
end
lonGrid = bestLon; latGrid = bestLat; residual = bestScore;
end

function grid = readFloatRaster(path, width, imageLength, byteOrder)
fid = fopen(path,'r',byteOrder);
if fid < 0
    error('PHASE_StaMPS:gridDiagnosticGeoRead','Cannot open %s.',path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
grid = fread(fid,[width imageLength],'single=>double')';
end

function [azStep, rgStep] = gridSteps(imageLength, width, maximum)
azStep = max(1,ceil(imageLength/maximum));
rgStep = max(1,ceil(width/maximum));
end

function extent = bufferedExtent(lon, lat, bufferMetres)
valid = isfinite(lon) & isfinite(lat);
lon = lon(valid); lat = lat(valid);
if isempty(lon)
    error('PHASE_StaMPS:gridDiagnosticCoordinates','No finite candidate coordinates exist.');
end
centreLat = mean([min(lat) max(lat)]);
latBuffer = bufferMetres / 111320;
lonScale = max(cosd(centreLat),0.01);
lonBuffer = bufferMetres / (111320*lonScale);
extent = [min(lon)-lonBuffer max(lon)+lonBuffer ...
    min(lat)-latBuffer max(lat)+latBuffer];
end

function createDaFigure(path, geo, candidates, extent, azStep, rgStep, threshold, runName, mode)
fig = figure('Visible','off','Color','w','Position',[50 50 1600 1000]);
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
ax = geoaxes(fig); hold(ax,'on'); setBasemap(ax);
drawRadarGrid(ax,geo,azStep,rgStep);
valid = isfinite(candidates.lon) & isfinite(candidates.lat) & isfinite(candidates.da);
geoscatter(ax,candidates.lat(valid),candidates.lon(valid),8,candidates.da(valid),'filled');
colorbar(ax); colormap(ax,parula(256));
geolimits(ax,extent(3:4),extent(1:2));
title(ax,{sprintf('Dispersione d''ampiezza dei candidati - run %s, D_A <= %.4g',runName,threshold), ...
    gridSubtitle(azStep,rgStep,mode)},'Interpreter','none');
exportgraphics(fig,path,'Resolution',220);
end

function createStageFigure(path, geo, candidates, extent, azStep, rgStep, threshold, runName, mode)
fig = figure('Visible','off','Color','w','Position',[50 50 1600 1000]);
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
ax = geoaxes(fig); hold(ax,'on'); setBasemap(ax);
drawRadarGrid(ax,geo,azStep,rgStep);
colours = [0.55 0.55 0.55; 0.95 0.50 0.08; 0.78 0.12 0.55];
names = {'Eliminato in ps_select','Eliminato al weeding', ...
    'Sopravvive al weeding (patch-local ps2)'};
handles = gobjects(0); labels = {};
for stage = 1:3
    use = candidates.stage == stage & isfinite(candidates.lon) & isfinite(candidates.lat);
    handles(end+1) = geoscatter(ax,candidates.lat(use),candidates.lon(use), ...
        9,colours(stage,:),'filled'); %#ok<AGROW>
    labels{end+1} = sprintf('%s: %d',names{stage},nnz(use)); %#ok<AGROW>
end
legend(ax,handles,labels,'Location','best','Interpreter','none');
geolimits(ax,extent(3:4),extent(1:2));
title(ax,{sprintf('Stadio massimo raggiunto dai candidati StaMPS - run %s, D_A <= %.4g',runName,threshold), ...
    gridSubtitle(azStep,rgStep,mode)},'Interpreter','none');
exportgraphics(fig,path,'Resolution',220);
end

function subtitle = gridSubtitle(azStep,rgStep,mode)
subtitle = sprintf(['Griglia di campionamento radar: linee ogni %d celle azimuth ', ...
    'e %d celle range; geolocalizzazione %s'],azStep,rgStep,strrep(mode,'_',' '));
end

function setBasemap(ax)
try
    geobasemap(ax,'satellite');
catch
    try, geobasemap(ax,'streets-light'); catch, end
end
end

function drawRadarGrid(ax, geo, azStep, rgStep)
azBoundaries = unique([-0.5:azStep:(geo.length-0.5), geo.length-0.5]);
rgBoundaries = unique([-0.5:rgStep:(geo.width-0.5), geo.width-0.5]);
azCurve = linspace(-0.5,geo.length-0.5,geo.length+1);
rgCurve = linspace(-0.5,geo.width-0.5,geo.width+1);
for k = 1:numel(azBoundaries)
    [lon,lat] = evaluateGeo(geo,azBoundaries(k)*ones(size(rgCurve)),rgCurve);
    geoplot(ax,lat,lon,'Color',[0.92 0.92 0.92],'LineWidth',0.35, ...
        'HandleVisibility','off');
end
for k = 1:numel(rgBoundaries)
    [lon,lat] = evaluateGeo(geo,azCurve,rgBoundaries(k)*ones(size(azCurve)));
    geoplot(ax,lat,lon,'Color',[0.92 0.92 0.92],'LineWidth',0.35, ...
        'HandleVisibility','off');
end
end

function [lon,lat] = evaluateGeo(geo,azBoundary,rgBoundary)
% Candidate indices are zero-based centres. Convert them to the one-based
% centre coordinate system used by MATLAB interpolation. Boundaries remain
% at +/- 0.5 around each integer-valued candidate centre.
azCentreCoordinate = azBoundary + 1;
rgCentreCoordinate = rgBoundary + 1;
if strcmp(geo.mode,'grid')
    lon = geo.lonF(azCentreCoordinate,rgCentreCoordinate);
    lat = geo.latF(azCentreCoordinate,rgCentreCoordinate);
else
    lon = geo.lon(rgBoundary,azBoundary);
    lat = geo.lat(rgBoundary,azBoundary);
end
end

function name = safeRunName(exportName, workDir)
name = regexprep(strtrim(char(string(exportName))), '[^A-Za-z0-9._-]+', '_');
if isempty(name)
    [~,name] = fileparts(workDir);
end
end

function writeJson(path, value)
fid = fopen(path,'w','n','UTF-8');
if fid < 0
    error('PHASE_StaMPS:gridDiagnosticMetadata','Cannot create %s.',path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
