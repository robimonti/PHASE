function throwIfStopped(engine)
%THROWIFSTOPPED Yield to the UI and stop at the next safe checkpoint.

drawnow limitrate
if ~isempty(engine) && isvalid(engine) && engine.StopRequested
    error('PHASE_Model_beta:hardStopped', ...
        'Processing was stopped by the user. Partial output remains available for diagnostics.');
end
end
