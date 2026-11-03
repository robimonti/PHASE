function info=covarianceReview(settings,meanInfo)
%COVARIANCEREVIEW Descriptive fit checks; never certify stationarity from a fit.
info=struct('status','diagnostics available','certifiesStationarity',false, ...
    'scope','Weighted discrepancy and residual-mean checks only. Conditional kernel fit, not a stationarity test.');
for axis={'Time','Space'}
    name=axis{1};lag=settings.(['lag' name]);empirical=settings.(['covariance' name]);
    fun=settings.([lower(name) 'Function']);curve=fun(settings.([lower(name) 'Parameters']),lag);
    count=settings.(['pairCount' name]);valid=lag>0 & isfinite(empirical) & count>0;
    scale=max(abs(empirical(1)),eps);
    info.([lower(name) 'RelativeRMSE'])=sqrt(sum(count(valid).*(empirical(valid)-curve(valid)).^2)/max(sum(count(valid)),eps))/scale;
end
info.meanToRMS=abs(meanInfo.mean)/max(meanInfo.rms,eps);
if info.meanToRMS>.1 || info.timeRelativeRMSE>.5 || info.spaceRelativeRMSE>.5
    info.status='review residual structure / covariance fit';
    warning('PHASE:residualCovarianceReview','Residual/covariance diagnostics need review: mean/RMS %.2f, temporal relative fit error %.2f, spatial %.2f. Detrending does not certify stationarity.', ...
        info.meanToRMS,info.timeRelativeRMSE,info.spaceRelativeRMSE);
end
end
