function exportSavedCovariance(metadataFile,outputFolder)
%EXPORTSAVEDCOVARIANCE Recover covariance plots without refitting a saved run.
metadata=jsondecode(fileread(metadataFile));
if ~isfield(metadata,'residualSettings')||~isfield(metadata.residualSettings,'timeModel')
    error('PHASE:customDiagnostics','Select custom_model.json from a native STC run.');
end
settings=metadata.residualSettings;
settings.timeFunction=kernel(settings.timeModel);
settings.spaceFunction=kernel(settings.spaceModel);
phase_model_custom.exportCovarianceDiagnostics(settings,outputFolder);
end
function evaluator=kernel(name)
switch name
    case 'gaussian',evaluator=@(p,d)p(1).*exp(-p(2).*d.^2);
    case 'exponential',evaluator=@(p,d)p(1).*exp(-p(2).*d);
    case {'gaussian+cosine','gaussian_cos','gaussian_cosine'}
        evaluator=@(p,d)p(1).*exp(-p(2).*d.^2).*cos(p(3).*d.*pi/2);
    otherwise,error('PHASE:customDiagnostics','Unknown saved covariance model: %s',name);
end
end
