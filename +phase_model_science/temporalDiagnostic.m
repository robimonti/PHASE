function info=temporalDiagnostic(times,residuals,folder,name)
%TEMPORALDIAGNOSTIC Inspect residual means, never fit a covariance curve as a trend.
% Residuals: PS x epochs. Block validation is a descriptive selection heuristic,
% not an independent-sample significance test or proof of stationarity.
times=double(times(:));means=mean(residuals,1,'omitnan')';
counts=sum(isfinite(residuals),1)';good=isfinite(means)&counts>0;
t=times(good);y=means(good);scale=max(range(t),1);origin=min(times);
u=(t-origin)/scale;score=inf(1,3);coefficients=cell(1,3);
if numel(t)>=8
    blocks=min(5,floor(numel(t)/3));labels=min(blocks,ceil((1:numel(t))'*blocks/numel(t)));
    for degree=0:2
        A=u.^(0:degree);coefficients{degree+1}=A\y;
        prediction=nan(size(y));
        for block=1:blocks
            train=labels~=block;test=~train;
            if nnz(train)<=degree+1,continue;end
            prediction(test)=A(test,:)*(A(train,:)\y(train));
        end
        score(degree+1)=mean((y-prediction).^2);
    end
end
degree=0;
if all(isfinite(score)) && score(2)<0.8*score(1),degree=1;end
if all(isfinite(score)) && score(3)<0.8*score(1) && score(3)<0.9*score(2),degree=2;end
rmsValue=sqrt(mean(residuals(:).^2,'omitnan'));
drift=0;
if degree>0
    curve=(u.^(0:degree))*coefficients{degree+1};drift=range(curve);
    if drift<max(0.5,0.1*rmsValue),degree=0;end
end
info=struct('suggestedDegree',degree,'blockedCV_MSE',score, ...
    'mean_mm',mean(residuals(:),'omitnan'),'rms_mm',rmsValue, ...
    'fittedDrift_mm',drift,'timeOrigin_days',origin,'timeScale_days',scale, ...
    'interpretation','Descriptive blocked validation of the AOI mean: >=20% improvement over constant; quadratic also >=10% over linear; drift >= max(0.5 mm, 0.1 residual RMS). Not a stationarity test. Local trends may cancel in the AOI mean.');
% Inspect individual PS as well: opposite slopes can hide in the AOI mean.
pointSlope=nan(size(residuals,1),1);pointCurvature=pointSlope;
pointScores=nan(size(residuals,1),3);pointDrift=pointScores;
for k=1:size(residuals,1)
    valid=isfinite(residuals(k,:));
    if nnz(valid)<4,continue;end
    localTime=(times(valid)-origin)/scale;
    polynomial=(localTime.^(0:2))\residuals(k,valid)';
    pointSlope(k)=polynomial(2)/scale*365.25;
    pointCurvature(k)=polynomial(3)/scale^2*365.25^2;
    if nnz(valid)<8,continue;end
    localY=residuals(k,valid)';
    blocks=min(5,floor(numel(localY)/3));
    labels=min(blocks,ceil((1:numel(localY))'*blocks/numel(localY)));
    for d=0:2
        A=localTime.^(0:d);prediction=nan(size(localY));
        for block=1:blocks
            train=labels~=block;test=~train;
            prediction(test)=A(test,:)*(A(train,:)\localY(train));
        end
        pointScores(k,d+1)=mean((localY-prediction).^2);
        pointDrift(k,d+1)=range(A*(A\localY));
    end
end
% Aggregate held-out errors across PS, rather than the AOI mean, so equal
% and opposite local drifts remain visible. Selection is descriptive only.
localScore=mean(pointScores,1,'omitnan');localDegree=0;
if all(isfinite(localScore)) && localScore(2)<.8*localScore(1),localDegree=1;end
if all(isfinite(localScore)) && localScore(3)<.8*localScore(1) && localScore(3)<.9*localScore(2),localDegree=2;end
localDrift=0;
if localDegree>0
    localDrift=sqrt(mean(pointDrift(:,localDegree+1).^2,'omitnan'));
    if localDrift<max(.5,.1*rmsValue),localDegree=0;end
end
info.suggestedPointDegree=localDegree;
info.pointBlockedCV_MSE=localScore;info.pointDriftRMS_mm=localDrift;
info.pointSlopePercentiles_mmYear=prctile(pointSlope,[10 50 90]);
info.pointCurvaturePercentiles_mmYear2=prctile(pointCurvature,[10 50 90]);
if nargin<3||isempty(folder),return;end
if ~isfolder(folder),mkdir(folder);end
if nargin<4,name='temporal_residual';end
linear=nan(size(times));quadratic=linear;
if ~isempty(coefficients{2}),linear=((times-origin)/scale).^(0:1)*coefficients{2};end
if ~isempty(coefficients{3}),quadratic=((times-origin)/scale).^(0:2)*coefficients{3};end
writetable(table(times,means,counts,linear,quadratic,'VariableNames', ...
    {'Time_days','ResidualMean_mm','PS_count','LinearFit_mm','QuadraticFit_mm'}),fullfile(folder,[name '.csv']));
writetable(table((1:size(residuals,1))',pointSlope,pointCurvature,'VariableNames', ...
    {'PS_index','LinearCoefficientAtStart_mmYear','QuadraticCoefficient_mmYear2'}), ...
    fullfile(folder,[name '_by_PS.csv']));
phase_project.writeJson(fullfile(folder,[name '.json']),info);
end
