function [retry,protectedDegree,polyDegree]=nativeTrendCheck(times,residuals,options,protectedDegree,polyDegree,folder,name)
%NATIVETRENDCHECK Reuse the automatic engine's polynomial and uncertainty path.
% At most one optional retry: add/protect temporal terms in the existing basis.
% This is model selection, not an edit to observations or the outlier mask.
info=phase_model_science.temporalDiagnostic(times,residuals,folder,name);
requested=0;
switch options.temporalTrend
    case 'auto',requested=max(info.suggestedDegree,info.suggestedPointDegree);
    case 'linear',requested=1;
    case 'quadratic',requested=2;
end
retry=protectedDegree==0 && requested>0;
if retry
    protectedDegree=requested;polyDegree=max(polyDegree,requested);
    if strcmp(options.temporalTrend,'auto') && info.suggestedPointDegree>0
        % Permit spatially varying slopes/curvatures using the native basis.
        polyDegree=max(polyDegree,min(3,requested+1));
    end
    fprintf('Residual temporal check: refitting native polynomial with degree %d temporal terms retained.\n',requested);
end
end
