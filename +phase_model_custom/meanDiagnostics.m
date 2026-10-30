function info=meanDiagnostics(c,r,ids,context,folder)
%MEANDIAGNOSTICS Descriptive checks, with no iid significance claim.
info=struct('mean',mean(r),'rms',sqrt(mean(r.^2)), ...
    'rawSecondMoment',mean(r.^2),'centeredSecondMoment',mean((r-mean(r)).^2), ...
    'interpretation','Finite-sample descriptive diagnostics; neither proof nor rejection of a zero-mean stochastic process.');
[times,~,groups]=unique(c(:,4));
means=accumarray(groups,r,[],@mean); counts=accumarray(groups,1);
writetable(table(times,means,counts,'VariableNames',{'Time_days','ResidualMean_mm','Count'}), ...
    fullfile(folder,'custom_residual_mean_by_epoch.csv'));
[point,~,groups]=unique(ids);means=accumarray(groups,r,[],@mean);
writetable(table(point,context.xy(point,1),context.xy(point,2),means, ...
    'VariableNames',{'PS_index','X_m','Y_m','ResidualMean_mm'}),fullfile(folder,'custom_residual_mean_by_PS.csv'));
info.warningRatio=0.1;
if abs(info.mean)>info.warningRatio*max(info.rms,eps)
    warning('PHASE:customResidualMean','Residual mean is %.4g mm (RMS %.4g mm). This descriptive warning does not trigger automatic detrending.',info.mean,info.rms);
end
phase_project.writeJson(fullfile(folder,'custom_residual_mean.json'),info);
end
