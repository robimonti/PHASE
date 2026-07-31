function [errors,warnings] = validateConfig(config,forRun)
%VALIDATECONFIG Validate visible Model settings before save or execution.

if nargin < 2, forRun = false; end
errors = {};
warnings = {};

if isempty(strtrim(config.filepathIN))
    if forRun, errors{end+1} = 'Select the displacement time-series XLSX/CSV file.'; end %#ok<AGROW>
elseif forRun && ~isfile(config.filepathIN)
    errors{end+1} = ['Input time-series file does not exist: ' config.filepathIN]; %#ok<AGROW>
end
if isempty(strtrim(config.pythonPath))
    errors{end+1} = 'Select a Python installation.'; %#ok<AGROW>
elseif forRun && ispc && ~isfile(config.pythonPath) && ...
        ~any(strcmpi(config.pythonPath,{'python','python3','py -3'}))
    errors{end+1} = ['Python executable does not exist: ' config.pythonPath]; %#ok<AGROW>
end
if ~isdatetime(config.t0IN) || isnat(config.t0IN)
    errors{end+1} = 'First displacement date is invalid.'; %#ok<AGROW>
end
if ~isfinite(config.markerSize) || config.markerSize <= 0
    errors{end+1} = 'Plot marker size must be greater than zero.'; %#ok<AGROW>
end

if config.flag_AOIbb
    polygon = config.aoi_polygon_lonlat;
    if ~isempty(polygon)
        if size(polygon,2) ~= 2 || size(polygon,1) < 3 || ...
                any(~isfinite(polygon(:))) || any(abs(polygon(:,1)) > 180) || ...
                any(abs(polygon(:,2)) > 90)
            errors{end+1} = 'The map AOI must contain at least three valid longitude/latitude vertices.'; %#ok<AGROW>
        else
            openPolygon = polygon;
            if size(openPolygon,1)>1 && isequal(openPolygon(1,:),openPolygon(end,:))
                openPolygon(end,:) = [];
            end
            if size(unique(openPolygon,'rows'),1) < 3 || ...
                    abs(polyarea(openPolygon(:,1),openPolygon(:,2))) < eps
                errors{end+1} = 'The map AOI must enclose a non-zero area.'; %#ok<AGROW>
            end
        end
    else
        bounds = [config.lonMinAOI config.lonMaxAOI config.latMinAOI config.latMaxAOI];
        if any(~isfinite(bounds))
            errors{end+1} = 'Enter all four geographic AOI bounds or draw a polygon on the map.'; %#ok<AGROW>
        elseif config.lonMinAOI >= config.lonMaxAOI || ...
                config.latMinAOI >= config.latMaxAOI || ...
                config.lonMinAOI < -180 || config.lonMaxAOI > 180 || ...
                config.latMinAOI < -90 || config.latMaxAOI > 90
            errors{end+1} = 'AOI bounds are invalid or outside longitude/latitude limits.'; %#ok<AGROW>
        end
    end
elseif isempty(strtrim(config.filepathAOI))
    if forRun, errors{end+1} = 'Select an AOI shapefile or enable the bounding box.'; end %#ok<AGROW>
elseif forRun && ~isfile(config.filepathAOI)
    errors{end+1} = ['AOI shapefile does not exist: ' config.filepathAOI]; %#ok<AGROW>
end

if strcmp(config.projDim,'1D')
    if ~strcmp(config.procType,'temporal') && ...
            (~isfinite(config.cline_resolution) || config.cline_resolution <= 0)
        errors{end+1} = 'Centerline resolution must be greater than zero for 1D interpolation.'; %#ok<AGROW>
    end
else
    if ~strcmp(config.procType,'temporal') && ...
            (~isfinite(config.grid_resolution) || config.grid_resolution <= 0)
        errors{end+1} = 'Grid resolution must be greater than zero for 2D interpolation.'; %#ok<AGROW>
    end
end
if ~ismember(config.procType,{'temporal','temporal&NNI'}) && ...
        (~isfinite(config.minMonths) || config.minMonths <= 0)
    errors{end+1} = 'Minimum spline-filtering period must be greater than zero.'; %#ok<AGROW>
end
if ~ismember(config.procType,{'temporal','temporal&NNI'}) && ...
        (~isfinite(config.step_t_ST) || config.step_t_ST <= 0)
    errors{end+1} = 'Temporal step must be greater than zero.'; %#ok<AGROW>
end

if strcmp(config.varNoise_method,'manual') && ...
        (~isfinite(config.varNoise_manual) || config.varNoise_manual < 0)
    errors{end+1} = 'Manual a-priori variance must be zero or greater.'; %#ok<AGROW>
elseif strcmp(config.varNoise_method,'coherence')
    if isempty(strtrim(config.coherence_dir))
        errors{end+1} = 'Select the coherence folder.'; %#ok<AGROW>
    elseif forRun && ~isfolder(config.coherence_dir)
        errors{end+1} = ['Coherence folder does not exist: ' config.coherence_dir]; %#ok<AGROW>
    end
    if ~isfinite(config.num_looks) || config.num_looks <= 0
        errors{end+1} = 'Number of looks must be greater than zero.'; %#ok<AGROW>
    end
end
if strcmp(config.num_spl_method,'manual') && ...
        (~isfinite(config.num_spl_manual) || config.num_spl_manual < 2)
    errors{end+1} = 'Manual temporal splines must be at least 2.'; %#ok<AGROW>
end
if strcmp(config.lambda_method,'manual') && ...
        (~isfinite(config.lambda_manual) || config.lambda_manual < 0)
    errors{end+1} = 'Manual temporal lambda must be zero or greater.'; %#ok<AGROW>
end
if strcmp(config.coll_proc,'prediction') && ...
        (~isfinite(config.coll_step_est) || config.coll_step_est <= 0)
    errors{end+1} = 'Collocation prediction step must be greater than zero.'; %#ok<AGROW>
end

if ismember(config.procType,{'temporal','temporal&NNI'})
    errors = [errors validateTemporalThreshold( ...
        config,'min_period_days_method','min_period_days', ...
        @(value) value > 0,'Minimum Fourier period must be greater than zero.')]; %#ok<AGROW>
    errors = [errors validateTemporalThreshold( ...
        config,'min_coll_snr_method','min_coll_snr', ...
        @(value) value > 0,'Minimum collocation SNR must be greater than zero.')]; %#ok<AGROW>
    errors = [errors validateTemporalThreshold( ...
        config,'min_coll_corr_samples_method','min_coll_corr_samples', ...
        @(value) value >= 1,'Minimum correlation length must be at least one sample.')]; %#ok<AGROW>
    errors = [errors validateTemporalThreshold( ...
        config,'spline_min_knot_intervals_method','spline_min_knot_intervals', ...
        @(value) value >= 2,'Minimum spline knot spacing must be at least two intervals.')]; %#ok<AGROW>
    errors = [errors validateTemporalThreshold( ...
        config,'spline_max_fraction_method','spline_max_fraction', ...
        @(value) value > 0 && value < 0.5, ...
        'Maximum spline fraction must be greater than zero and below 0.5.')]; %#ok<AGROW>
end

if config.flag_tsExtr
    if isempty(strtrim(config.filepath_EXTR))
        errors{end+1} = 'Select the query-points file.'; %#ok<AGROW>
    elseif forRun && ~isfile(config.filepath_EXTR)
        errors{end+1} = ['Query-points file does not exist: ' config.filepath_EXTR]; %#ok<AGROW>
    end
end

if strcmp(config.procType,'spatialDET')
    if strcmp(config.projDim,'1D')
        errors = [errors validateDeterministic(config,'1D')]; %#ok<AGROW>
    else
        errors = [errors validateDeterministic(config,'2D')]; %#ok<AGROW>
    end
end
if strcmp(config.procType,'spatialSTC')
    if strcmp(config.projDim,'1D')
        pairs = {'dtCov_method_STC1D','dtCov_STC1D'; ...
            'dsCov_method_STC1D','dsCov_STC1D'};
    else
        pairs = {'dtCov_method_STC2D','dtCov_STC2D'; ...
            'dsCov_method_STC2D','dsCov_STC2D'};
    end
    for k = 1:size(pairs,1)
        if strcmp(config.(pairs{k,1}),'manual') && ...
                (~isfinite(config.(pairs{k,2})) || config.(pairs{k,2}) <= 0)
            errors{end+1} = 'Manual covariance bins must be greater than zero.'; %#ok<AGROW>
        end
    end
end

if ~forRun && isempty(errors) && isempty(strtrim(config.filepathIN))
    warnings{end+1} = 'The configuration can be saved now; an input file is required only before Start.'; %#ok<AGROW>
end
end

function errors = validateTemporalThreshold(config,methodName,valueName,predicate,message)
errors = {};
method = config.(methodName);
if ~any(strcmp(method,{'auto','manual'}))
    errors{end+1} = [strrep(methodName,'_',' ') ' must be Automatic or Manual.'];
elseif strcmp(method,'manual')
    value = config.(valueName);
    if ~isfinite(value) || ~predicate(value)
        errors{end+1} = message;
    end
end
end

function errors = validateDeterministic(config,dimension)
errors = {};
suffix = ['DET' dimension];
noiseMethod = config.(['varNoise_' suffix]);
if strcmp(noiseMethod,'manual')
    value = config.(['varNoise_manual_' suffix]);
    if ~isfinite(value) || value < 0
        errors{end+1} = [dimension ' manual deterministic noise must be zero or greater.']; %#ok<AGROW>
    end
end
splinesMethod = config.(['num_spl_method_' suffix]);
if strcmp(splinesMethod,'manual')
    values = [config.(['num_spl_row_manual_' suffix]) ...
        config.(['num_spl_col_manual_' suffix])];
    if strcmp(dimension,'2D')
        values(end+1) = config.num_spl_t_manual_DET2D;
    end
    if any(~isfinite(values) | values < 2)
        errors{end+1} = [dimension ' manual spline counts must all be at least 2.']; %#ok<AGROW>
    end
end
lambdaMethod = config.(['lambda_method_' suffix]);
if strcmp(lambdaMethod,'manual')
    value = config.(['lambda_manual_' suffix]);
    if ~isfinite(value) || value < 0
        errors{end+1} = [dimension ' manual lambda must be zero or greater.']; %#ok<AGROW>
    end
end
end
