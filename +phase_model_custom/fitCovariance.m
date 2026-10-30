function settings=fitCovariance(c,r,config,context,pointIds)
%FITCOVARIANCE Use the original empirical product-covariance functions.
[~,it]=ismember(c(:,4),context.times);
matrix=nan(size(context.xy,1),numel(context.times));
matrix(sub2ind(size(matrix),pointIds,it))=r;
if strcmp(config.projDim,'1D')
    settings=phase_model_custom.nativeCovariance1D(matrix,[context.s,zeros(size(context.s))],context.times,config);
else
    settings=phase_model_custom.nativeCovariance(matrix,context.xy,context.times,config);
end
settings.referenceDay=context.referenceDay;
settings.timeModel=config.(['tCovModel_STC' config.projDim]);
settings.spaceModel=config.(['sCovModel_STC' config.projDim]);
end
