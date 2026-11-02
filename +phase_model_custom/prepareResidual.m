function residual=prepareResidual(c,y,ids,fit,config)
%PREPARERESIDUAL Original residual algorithms, after one custom mean fit.
r=y-fit.mean;
if strcmp(config.procType,'spatialDET')
    residual=phase_model_custom.nativeDetFit(c,r,ids,fit,config); return
end
% Conservative working estimate: covariance, noise, Cholesky and temporaries.
estimatedGB=(6*8*numel(y)^2+128*2^20)/2^30;
if estimatedGB>config.customMemoryGB
    error('PHASE:customMemory','Dense STC needs approximately %.2f GiB for %d retained space-time samples; budget: %.2f GiB. Reduce the AOI or increase the budget if RAM permits.',estimatedGB,numel(y),config.customMemoryGB);
end
phase_model_custom.reportProgress(config,56,sprintf('Fitting empirical residual covariances — %d samples, estimated matrix RAM %.2f GiB',numel(y),estimatedGB));
settings=phase_model_custom.fitCovariance(c,r,config,fit.context,ids);
phase_model_custom.reportProgress(config,57,'Building the observation covariance matrix');
if strcmp(config.projDim,'1D')
    [C,~,N]=compute_C_obs_obs_nugget([c(:,4),c(:,3)],settings.timeFunction, ...
        settings.timeParameters,settings.spaceFunction,settings.spaceParameters,settings.noiseVariance,fit.context.referenceDay);
else
    [C,~,N]=compute_C_obs_obs_nugget_2D([c(:,4),c(:,1:2)],settings.timeFunction, ...
        settings.timeParameters,settings.spaceFunction,settings.spaceParameters,settings.noiseVariance,fit.context.referenceDay);
end
C=double(C); C=(C+C')/2;
phase_model_custom.reportProgress(config,58,'Factorising the covariance matrix for collocation');
[factor,bad]=chol(C,'lower'); jitter=0;
if bad
    jitter=1e-8*max(max(diag(C)),eps);
    for attempt=1:8
        [factor,bad]=chol(C+jitter*eye(size(C)),'lower');
        if ~bad, C=C+jitter*eye(size(C)); break; end
        jitter=jitter*10;
    end
end
if bad, error('PHASE:customCovariance','The native covariance fit is not positive definite. Review covariance settings.'); end
settings.numericalJitter=jitter;
residual=struct('kind','STC','native',true,'settings',settings,'C',C, ...
    'noiseMatrix',N,'factor',factor,'alpha',factor'\(factor\r));
end
