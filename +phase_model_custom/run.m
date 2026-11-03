function [lonlat,dates,times,atPS,onGrid,stdGrid,atQuery] = run( ...
    clean,raw,xy,datesIN,timesIN,centerline,xGrid,yGrid,xyAOI,xyQuery,utmZone,config,outputDir,filterFit)
%RUN Isolated custom workflow after the original unmodified outlier filter.
% 2) Prepare variables; 3) deterministic fit; 4) residual fit;
% 5) joint uncertainty propagation; 6) outputs and provenance.
fprintf('Custom model: two stages; automatic model reduction is bypassed.\n');
outputDir=fullfile(outputDir,'files','customModel');
if ~isfolder(outputDir), mkdir(outputDir); end
timesIN=double(timesIN(:)'); timesIN=timesIN-timesIN(1);
datesIN=datesIN(:);
if numel(timesIN)<3 || any(diff(timesIN)<=0)
    error('PHASE:customTime','Custom models require at least three strictly increasing acquisition dates.');
end
definition=phase_model_custom.parseDefinition(config);
phase_model_custom.validateOptions(config);
master=find(all(abs(raw)<1e-6,1));
if isempty(master)
    error('PHASE:customReference','No zero-displacement reference acquisition was found. Supply series with a common explicit reference epoch.');
end
if numel(master)>1
    [~,closest]=min(abs(timesIN(master)-mean(timesIN)));
    master=master(closest);
end
context=struct('referenceDay',timesIN(master),'origin',mean(xy,1));
context.times=timesIN; context.xy=xy; context.outputDir=outputDir;
context.utmZone=utmZone;context.coordinateSystem=config.customCoordinates;
if strcmp(config.customCoordinates,'lonlat')
    [lat,lon]=utm2deg(xy(:,1),xy(:,2),repmat(utmZone,size(xy,1),1));
    context.origin=mean([lon,lat],1);
end
if isfinite(config.customOriginX),context.origin(1)=config.customOriginX;end
if isfinite(config.customOriginY),context.origin(2)=config.customOriginY;end
context.filterResidual=raw-filterFit';
if isfield(config,'stopCheck'), context.stopCheck=config.stopCheck; end
np=size(xy,1); nt=numel(timesIN);
if strcmp(config.projDim,'1D')
    s=centerline.projected_distances(:);
    gridXY=centerline.xy_centerline;
    gridS=centerline.s(:);
else
    s=zeros(np,1); gridXY=[xGrid(:),yGrid(:)];
    inside=inpolygon(gridXY(:,1),gridXY(:,2),xyAOI(:,1),xyAOI(:,2));
    gridXY=gridXY(inside,:); gridS=zeros(size(gridXY,1),1);
end
if isempty(gridXY), error('PHASE:customGrid','No prediction points lie inside the AOI.'); end
allCoordinates=expand(xy,s,timesIN);
context.s=s;
context.domainCoordinates=[allCoordinates;expand(gridXY,gridS,[0,timesIN(end)])];
queryS=[];
if ~isempty(xyQuery)
    if strcmp(config.projDim,'1D'), queryS=projectQueries(xyQuery,gridXY,gridS);
    else, queryS=zeros(size(xyQuery,1),1); end
    context.domainCoordinates=[context.domainCoordinates;expand(xyQuery,queryS,[0,timesIN(end)])];
end
allIds=repmat((1:np)',nt,1);
valid=isfinite(clean(:));
c=allCoordinates(valid,:); y=clean(valid); ids=allIds(valid);
if numel(y)<max(5,2*nnz([definition.parameters.estimate])+1)
    error('PHASE:customObservations','Too few valid observations remain after outlier rejection.');
end
fprintf('Custom deterministic fit: %d observations, %d estimated parameters.\n',numel(y),nnz([definition.parameters.estimate]));
phase_model_custom.checkStop(config);
phase_model_custom.reportProgress(config,54,'Fitting the custom deterministic model');
fit=phase_model_custom.fitDeterministic(definition,c,y,context);
fit=phase_model_custom.augmentTemporalTrend(fit,c,y,ids,config,outputDir);
if ~isfolder(outputDir), mkdir(outputDir); end
customMean=phase_model_custom.evaluate(fit.definition,fit.values,c,fit.context);
residualForCovariance=y-fit.mean;
writetable(table(ids,c(:,4),y,customMean,fit.mean-customMean,residualForCovariance, ...
    'VariableNames',{'PS_index','Time_days','Observed_mm','CustomModel_mm','TemporalTrend_mm','ResidualUsed_mm'}), ...
    fullfile(outputDir,'custom_observation_components.csv'));
meanDiagnostics=phase_model_custom.meanDiagnostics(c,y-fit.mean,ids,context,outputDir);
fprintf('Custom residual fit: %s.\n',config.procType);
phase_model_custom.reportProgress(config,56,'Estimating residual structure and uncertainty');
residual=phase_model_custom.prepareResidual(c,y,ids,fit,config);
covarianceReview=struct();
if strcmp(residual.kind,'STC')
    covarianceReview=phase_model_science.covarianceReview(residual.settings,meanDiagnostics);
    phase_project.writeJson(fullfile(outputDir,'covariance_review.json'),covarianceReview);
end
if strcmp(residual.kind,'STC')
    residual.parameterCovariance=fit.L*residual.C*fit.L';
    residual.whitenedJacobian=residual.factor\fit.J;
end
if any(definition.assignedCovariance(:))
    [~,residual.assignedJacobian]=phase_model_custom.evaluate(fit.definition,fit.values,c,fit.context);
    if strcmp(residual.kind,'STC')
        residual.whitenedAssignedJacobian=residual.factor\(residual.assignedJacobian+fit.J*fit.assignedSensitivity);
    end
end
times=(0:config.step_t_ST:timesIN(end))';
dates=datesIN(1)+days(times);
q=expand(gridXY,gridS,times);
prediction=predictBatches(q,0,'Output grid',59,65);
onGrid=reshape(prediction.total,size(gridXY,1),[]);
stdGrid=reshape(sqrt(prediction.variance),size(gridXY,1),[]);
obsPrediction=predictBatches(allCoordinates,context.referenceDay,'Input PS',65,68);
atPS=reshape(obsPrediction.total,np,nt);
atQuery=[];queryPrediction=struct();queryCoordinates=zeros(0,4);
if ~isempty(xyQuery)
    if strcmp(config.projDim,'1D')
        queryS=projectQueries(xyQuery,gridXY,gridS);
    else
        queryS=zeros(size(xyQuery,1),1);
    end
    queryCoordinates=expand(xyQuery,queryS,times);
    queryPrediction=predictBatches(queryCoordinates,0,'Query points',68,69);
    atQuery=reshape(queryPrediction.total,size(xyQuery,1),[]);
    count=size(xyQuery,1);
    identifiers="q"+string((1:count)');
    if isfield(config,'queryPointIds') && numel(config.queryPointIds)==count
        identifiers=string(config.queryPointIds(:));
    end
    [latitude,longitude]=utm2deg(xyQuery(:,1),xyQuery(:,2),repmat(utmZone,count,1));
    queryDates=repelem(dates,count);queryDates.Format='yyyy-MM-dd';
    queryTable=table(repmat(identifiers,numel(times),1),queryDates, ...
        repmat(longitude,numel(times),1),repmat(latitude,numel(times),1),queryCoordinates(:,4), ...
        queryPrediction.deterministic,queryPrediction.residual,queryPrediction.total,sqrt(queryPrediction.variance), ...
        'VariableNames',{'ID','Date','Longitude','Latitude','Time_days','Deterministic_mm','Residual_mm','Total_mm','Std_mm'});
    queryTable.CustomModel_mm=queryPrediction.customDeterministic;
    queryTable.TemporalTrend_mm=queryPrediction.temporalTrend;
    writetable(queryTable,fullfile(outputDir,'custom_query_timeseries.csv'));
end
[lat,lon]=utm2deg(gridXY(:,1),gridXY(:,2),repmat(utmZone,size(gridXY,1),1));
lonlat=[lon,lat];
% Full parameter covariance includes correlations induced by both stages.
if strcmp(residual.kind,'STC')
    parameterCovariance=fit.L*residual.C*fit.L';
else
    parameterCovariance=residual.noiseVariance*(fit.L*fit.L');
end
allParameterCovariance=zeros(numel(fit.values));
nFree=numel(fit.free);
allParameterCovariance(fit.free,fit.free)=parameterCovariance(1:nFree,1:nFree);
assignedInfluence=eye(numel(fit.values));
assignedInfluence(fit.free,:)=fit.assignedSensitivity(1:nFree,:);
allParameterCovariance=allParameterCovariance+assignedInfluence*definition.assignedCovariance*assignedInfluence';
metadata=struct('workflow','custom two-stage','family',config.procType, ...
    'dimension',config.projDim,'definition',definition,'parameterValues',fit.values, ...
    'parameterCovariance',allParameterCovariance,'estimatedParameterIndices',fit.free, ...
    'originUTM',mean(xy,1),'timeOrigin',char(string(datesIN(1))), ...
    'inputReferenceDate',char(string(datesIN(master))), ...
    'outputReferenceDate',char(string(datesIN(1))), ...
    'units','LOS displacement mm; metric x,y,s; time days', ...
    'uncertainty','Local sampling/prediction uncertainty, conditional on selected covariance/spline settings and outlier mask; includes deterministic-residual cross covariance. Excludes covariance hyperparameter/mean-field selection uncertainty, regularisation bias and model discrepancy.');
metadata.temporalTrend=fit.temporalTrend;
metadata.covarianceReview=covarianceReview;
metadata.augmentedParameterNames=[{definition.parameters(fit.free).name}, ...
    cellstr("residualMeanCoefficient"+string(1:(size(fit.J,2)-numel(fit.free))))];
metadata.augmentedParameterCovariance=parameterCovariance+fit.assignedSensitivity*definition.assignedCovariance*fit.assignedSensitivity';
metadata.residualMeanPolicy=config.customResidualMean;metadata.meanDiagnostics=meanDiagnostics;
metadata.coordinateSystem=config.customCoordinates;metadata.modelCoordinateOrigin=context.origin;
metadata.modelTimeOrigin=char(string(datesIN(1)));
metadata.assignedParameterAssumption='Assigned parameter errors are independent of observation errors; default covariance zero.';
if strcmp(residual.kind,'STC')
    metadata.residualSettings=rmfield(residual.settings,{'timeFunction','spaceFunction'});
else
    metadata.splineKnotsX=residual.nativeState.knots_x;metadata.splineKnotsY=residual.nativeState.knots_y;
    metadata.splineLambda=residual.lambda; metadata.noiseVariance=residual.noiseVariance;
end
if ~isfolder(outputDir), mkdir(outputDir); end
save(fullfile(outputDir,'custom_model_results.mat'),'metadata','prediction', ...
    'obsPrediction','clean','raw','filterFit','gridXY','gridS','times','dates', ...
    'xy','timesIN','datesIN','valid','ids','xyQuery','queryCoordinates','queryPrediction','residualForCovariance','-v7.3');
phase_project.writeJson(fullfile(outputDir,'custom_model.json'),metadata);
copyfileIfNeeded(definition.functionFile,outputDir);
% Tidy, explicit component exports; the existing workbook is generated by
% ProcessingEngine in exactly its established format.
componentTable=table(q(:,1),q(:,2),q(:,3),q(:,4),prediction.deterministic, ...
    prediction.residual,prediction.total,sqrt(prediction.variance), ...
    prediction.deterministicVariance,prediction.residualVariance,prediction.crossCovariance, ...
    'VariableNames',{'X_m','Y_m','s_m','Time_days','Deterministic_mm', ...
    'Residual_mm','Total_mm','Std_mm','DeterministicVariance_mm2', ...
    'ResidualErrorVariance_mm2','CrossCovariance_mm2'});
componentTable.CustomModel_mm=prediction.customDeterministic;
componentTable.TemporalTrend_mm=prediction.temporalTrend;
writetable(componentTable,fullfile(outputDir,'custom_model_components.csv'));
fprintf('Custom model completed. Diagnostics: %s\n',fullfile(outputDir,'custom_model.json'));
    function p=predictBatches(query,reference,label,startPercent,endPercent)
        % Bound temporary matrices while allowing efficient BLAS-sized blocks.
        p=struct(); batch=max(1,min(1024,floor(128*2^20/(8*numel(y)*8))));
        clock=tic;
        phase_model_custom.reportProgress(config,startPercent,sprintf('%s: residual prediction and uncertainty (0 / %d)',label,size(query,1)));
        for first=1:batch:size(query,1)
            last=min(first+batch-1,size(query,1));
            part=phase_model_custom.predict(query(first:last,:),reference,c,ids,fit,residual);
            names=fieldnames(part);
            for j=1:numel(names)
                name=names{j};
                if first==1, p.(name)=zeros(size(query,1),1); end
                p.(name)(first:last)=part.(name);
            end
            drawnow limitrate
            phase_model_custom.checkStop(config);
            if toc(clock)>=2 || last==size(query,1)
                fraction=last/size(query,1);
                phase_model_custom.reportProgress(config,startPercent+(endPercent-startPercent)*fraction, ...
                    sprintf('%s: residual prediction and uncertainty — %.0f%% (%d / %d space-time samples)',label,100*fraction,last,size(query,1)));
                clock=tic;
            end
        end
    end
end
function c=expand(xy,s,t)
n=size(xy,1); c=[repmat(xy,numel(t),1),repmat(s(:),numel(t),1),kron(t(:),ones(n,1))];
end
function s=projectQueries(q,xy,chainage)
s=zeros(size(q,1),1);
for k=1:size(q,1)
    a=xy(1:end-1,:); delta=diff(xy);
    fraction=max(0,min(1,sum((q(k,:)-a).*delta,2)./max(sum(delta.^2,2),eps)));
    projected=a+fraction.*delta;
    [~,j]=min(sum((projected-q(k,:)).^2,2));
    s(k)=chainage(j)+fraction(j)*(chainage(j+1)-chainage(j));
end
end
function copyfileIfNeeded(file,outputDir)
if isempty(file), return; end
[~,name,ext]=fileparts(file);
folder=fullfile(outputDir,'custom_model_source'); if ~isfolder(folder), mkdir(folder); end
destination=fullfile(folder,[name ext]);
if ~strcmp(file,destination), copyfile(file,destination); end
end
