function stop=checkStop(options)
%CHECKSTOP Poll the host cancellation request, including inside optimisers.
stop=false;
if isfield(options,'stopCheck') && ~isempty(options.stopCheck)
    options.stopCheck();
end
end
