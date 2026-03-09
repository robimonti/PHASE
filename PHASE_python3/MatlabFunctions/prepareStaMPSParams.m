function paramsData = prepareStaMPSParams(paramsStruct)
    
% prepareStaMPSParams Extracts StaMPS parameters from a struct and formats them for the report


paramsData = {};
paramNames = fieldnames(paramsStruct); % Get all variable names in the struct

% Find the index of 'heading' and include only variables from there onward
startIdx = find(strcmp(paramNames, 'heading'), 1); 

% If 'heading' is found, process only from 'heading' onwards
if ~isempty(startIdx)
    paramNames = paramNames(startIdx:end); % Filter paramNames to start from 'heading'
else
    % If 'heading' is not found, return an empty cell array
    paramsData = {};
    return;
end

% Loop through the filtered list of parameters
for i = 1:length(paramNames)
    paramName = paramNames{i};
    paramValue = paramsStruct.(paramName); % Get value of the parameter

    % Convert parameter to a string for readability
    if isnumeric(paramValue)
        paramValueStr = num2str(paramValue);
    elseif ischar(paramValue)
        paramValueStr = paramValue;
    else
        paramValueStr = 'Unsupported Type';
    end

    % Add name-value pair to the cell array
    paramsData = [paramsData; {sprintf('%s: %s', paramName, paramValueStr), ''}];
end
end