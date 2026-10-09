function [timeDays, displacement] = insertMasterEpoch(days, masterDay, firstDay, series)
%INSERTMASTEREPOCH Add the zero-displacement master on a valid date axis.
days = double(days(:)');
if size(series,2) ~= numel(days) || any(~isfinite(days)) || ...
        ~isscalar(masterDay) || ~isfinite(masterDay)
    error('PHASE_StaMPS:timeAxisInvalid', ...
        'StaMPS time-series columns and acquisition dates do not match.');
end
if any(diff(days) <= 0)
    error('PHASE_StaMPS:timeAxisOrder', ...
        'StaMPS acquisition dates must be strictly increasing.');
end
if any(days == masterDay)
    timeDays = days - firstDay;
    displacement = series;
    return
end
before = nnz(days < masterDay);
timeDays = [days(1:before),masterDay,days(before+1:end)] - firstDay;
displacement = [series(:,1:before),zeros(size(series,1),1), ...
    series(:,before+1:end)];
end
