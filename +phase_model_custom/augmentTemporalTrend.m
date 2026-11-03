function fit=augmentTemporalTrend(fit,c,y,ids,config,folder)
%AUGMENTTEMPORALTREND A second, optional mean fit with shared-data derivatives.
% beta = M*(y-m(theta(y))); D(beta)/D(y) = M - M*J*L.
% Appending these derivatives preserves cross covariance with both the custom
% fit and the residual prediction; fixed custom parameters remain fixed.
options=phase_model_science.options(config);
[~,it]=ismember(c(:,4),fit.context.times);
matrix=nan(size(fit.context.xy,1),numel(fit.context.times));
matrix(sub2ind(size(matrix),ids,it))=y-fit.mean;
diagnostic=phase_model_science.temporalDiagnostic(fit.context.times,matrix,folder,'temporal_before');
automatic=~isfield(config,'customResidualField') || strcmp(config.customResidualField,'automatic');
if automatic
    [trend,B,M]=phase_model_science.selectResidualField(c,y-fit.mean,ids,config.projDim,fit.context.referenceDay,config);
    trend.diagnostic=diagnostic;trend.scale=trend.timeScale;
    fit.temporalTrend=trend;
    degree=trend.degree;
    if degree==0 && diagnostic.suggestedDegree>0
        % Retain the independently diagnosed common mean if a spatial field
        % cannot generalise. Do not replace it with an unsupported local fit.
        trend.degree=diagnostic.suggestedDegree;trend.spatialDegree=0;
        trend.spatialKind='polynomial';trend.spatialIntervals=0;trend.lambda=0;
        trend.reason=['Common-mean fallback. ' trend.reason];
        B=phase_model_science.spaceTimeBasis(c,trend,fit.context.referenceDay);
        M=pinv(B);trend.coefficients=M*(y-fit.mean);
        fit.temporalTrend=trend;degree=trend.degree;
    end
    if degree==0
        fprintf('Residual field: %s\n',trend.reason);
        return
    end
    beta=trend.coefficients;
    if isfield(trend,'spatialPredictionNeedsReview') && trend.spatialPredictionNeedsReview
        warning('PHASE:residualFieldCoverage','Observed residual drift was corrected, but prediction in withheld spatial regions needs review. Output uncertainty is conditional on the selected mean field; spatial extrapolation/model bias is not included.');
    end
else
    degree=0;
    switch options.temporalTrend
        case 'auto',degree=diagnostic.suggestedDegree;
        case 'linear',degree=1;
        case 'quadratic',degree=2;
    end
    fit.temporalTrend=struct('degree',degree,'policy',options.temporalTrend, ...
        'diagnostic',diagnostic,'scale',max(range(fit.context.times),1),'coefficients',[]);
    if degree==0,return;end
    B=phase_model_science.trendBasis(c(:,4),fit.context.referenceDay,fit.temporalTrend.scale,degree);
    if rank(B)<degree,error('PHASE:temporalTrend','Insufficient epochs for the residual polynomial.');end
    M=pinv(B);
    beta=M*(y-fit.mean);
end
Ltrend=M-(M*fit.J)*fit.L;
[~,assignedJ]=phase_model_custom.evaluate(fit.definition,fit.values,c,fit.context);
assignedTrend=-M*(assignedJ+fit.J*fit.assignedSensitivity);
fit.mean=fit.mean+B*beta;
fit.J=[fit.J B];fit.L=[fit.L;Ltrend];
fit.assignedSensitivity=[fit.assignedSensitivity;assignedTrend];
fit.temporalTrend.coefficients=beta;
matrix(sub2ind(size(matrix),ids,it))=y-fit.mean;
fit.temporalTrend.after=phase_model_science.temporalDiagnostic(fit.context.times,matrix,folder,'temporal_after');
fprintf('Residual deterministic correction: temporal degree %d, %d coefficients; subtracted before residual fitting and restored in predictions.\n',degree,numel(beta));
end
