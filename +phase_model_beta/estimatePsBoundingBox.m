function bounds = estimatePsBoundingBox(filepath)
%ESTIMATEPSBOUNDINGBOX Read the complete PS extent from PHASE XLSX/CSV.

filepath = char(string(filepath));
if isempty(strtrim(filepath)) || ~isfile(filepath)
    error('PHASE_Model_beta:inputMissing', ...
        'Select an existing displacement time-series XLSX/CSV file first.');
end
[~,~,extension] = fileparts(filepath);
if ~any(strcmpi(extension,{'.xlsx','.csv'}))
    error('PHASE_Model_beta:invalidInputType', ...
        'The displacement time series must be an XLSX or CSV file.');
end

if strcmpi(extension,'.csv')
    values = readmatrix(filepath,'DecimalSeparator',',');
else
    values = readmatrix(filepath);
end
if size(values,1) < 3 || size(values,2) < 3
    error('PHASE_Model_beta:invalidInputLayout', ...
        'The time-series file does not contain PHASE longitude/latitude columns.');
end

coordinates = values(3:end,2:3);
valid = all(isfinite(coordinates),2) & ...
    coordinates(:,1) >= -180 & coordinates(:,1) <= 180 & ...
    coordinates(:,2) >= -90 & coordinates(:,2) <= 90;
coordinates = coordinates(valid,:);
if isempty(coordinates)
    error('PHASE_Model_beta:noCoordinates', ...
        'No valid PS longitude/latitude coordinates were found.');
end

lonLimits = [min(coordinates(:,1)) max(coordinates(:,1))];
latLimits = [min(coordinates(:,2)) max(coordinates(:,2))];
lonPadding = max(diff(lonLimits) * 0.005,1e-6);
latPadding = max(diff(latLimits) * 0.005,1e-6);
bounds = [ ...
    max(-180,lonLimits(1)-lonPadding), ...
    min(180,lonLimits(2)+lonPadding), ...
    max(-90,latLimits(1)-latPadding), ...
    min(90,latLimits(2)+latPadding)];
end
