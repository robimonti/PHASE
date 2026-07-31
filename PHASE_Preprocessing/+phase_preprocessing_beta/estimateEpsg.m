function code = estimateEpsg(lonMin, latMin, lonMax, latMax)
%ESTIMATEEPSG Select a projected WGS 84 CRS from the AOI centre.
%
% UTM is used between 80 S and 84 N, including the standard Norway and
% Svalbard zone exceptions. Polar stereographic CRSs are used beyond UTM.

values = double([lonMin latMin lonMax latMax]);
if any(~isfinite(values)) || lonMin >= lonMax || latMin >= latMax || ...
        lonMin < -180 || lonMax > 180 || latMin < -90 || latMax > 90
    error('PHASE_Preprocessing_beta:invalidAOIForCRS', ...
        'A valid AOI is required to estimate the output CRS.');
end

longitude = (lonMin + lonMax) / 2;
latitude = (latMin + latMax) / 2;
if latitude >= 84
    code = 3413; % WGS 84 / NSIDC Sea Ice Polar Stereographic North
    return
elseif latitude <= -80
    code = 3031; % WGS 84 / Antarctic Polar Stereographic
    return
end

zone = min(60,max(1,floor((longitude + 180) / 6) + 1));
if latitude >= 56 && latitude < 64 && longitude >= 3 && longitude < 12
    zone = 32;
elseif latitude >= 72 && latitude < 84
    if longitude >= 0 && longitude < 9
        zone = 31;
    elseif longitude < 21
        zone = 33;
    elseif longitude < 33
        zone = 35;
    elseif longitude < 42
        zone = 37;
    end
end

if latitude >= 0
    code = 32600 + zone;
else
    code = 32700 + zone;
end
end
