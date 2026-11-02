function reportProgress(config,percentage,message)
%REPORTPROGRESS Send custom stages to the host monitor; log for script users.
if isfield(config,'progressCallback') && ~isempty(config.progressCallback)
    config.progressCallback(percentage,message);
else
    fprintf('%s\n',message);
end
end
