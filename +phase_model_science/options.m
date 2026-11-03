function options=options(config)
%OPTIONS Shared, opt-in residual controls for automatic/custom DET and STC.
options=struct('temporalTrend','diagnose','correlationSpace',NaN, ...
    'correlationTime',NaN,'correlationTolerance',25, ...
    'correlationSpaceMax',NaN,'correlationTimeMax',NaN);
names=fieldnames(options);
for k=1:numel(names)
    if isfield(config,names{k}),options.(names{k})=config.(names{k});end
end
if ~ismember(options.temporalTrend,{'diagnose','auto','linear','quadratic'})
    error('PHASE:residualOptions','Select diagnose, auto, linear or quadratic temporal trend.');
end
for name={'correlationSpace','correlationTime','correlationSpaceMax','correlationTimeMax'}
    value=options.(name{1});
    if ~isnumeric(value)||~isscalar(value)||(~isnan(value)&&(~isfinite(value)||value<=0))
        error('PHASE:correlationRange','Correlation lengths must be positive or empty for the native fit.');
    end
end
v=options.correlationTolerance;
if ~isnumeric(v)||~isscalar(v)||~isfinite(v)||v<=0||v>=100
    error('PHASE:correlationRange','Correlation-length refinement must be between 0 and 100 percent.');
end
for dimension={'Space','Time'}
    prior=options.(['correlation' dimension{1}]);limit=options.(['correlation' dimension{1} 'Max']);
    if isfinite(prior)&&isfinite(limit)&&limit<prior
        error('PHASE:correlationRange','Maximum correlation length must not be below its suggested value.');
    end
end
end
