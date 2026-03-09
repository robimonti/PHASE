function setupPythonEnvironment_app(pythonPath_input)

% setupPythonEnvironment Automated Python Environment Setup for MATLAB

% Check the current Python environment
pyInfo = pyenv();

% If the environment is already loaded and functional
if pyInfo.Status == "Loaded"
    disp('Using the current Python environment:');
    disp(pyInfo.Executable);
else
    % Prompt the user to provide the Python path
    fprintf('MATLAB is not currently using a valid Python environment.\n');
    pythonPath = pythonPath_input;
    
    % Set the provided Python environment
    try
        pyenv('Version', pythonPath);
        fprintf('Successfully configured Python environment: %s\n', pythonPath);
    catch ME
        error('Failed to configure Python environment: %s\n', ME.message);
    end
end

% Validate the required Python package
try
    py.importlib.import_module('openpyxl');
    disp('The required package "openpyxl" is available.');
catch
    % Install the package if not found
    fprintf('The "openpyxl" package is not installed in the current Python environment.\n');
    choice = input('Would you like to install it automatically? (y/n): ', 's');
    if lower(choice) == 'y'
        try
            system([pyenv().Executable, ' -m pip install openpyxl']);
            disp('Successfully installed "openpyxl".');
        catch ME
            error('Failed to install "openpyxl": %s\n', ME.message);
        end
    else
        error('Cannot proceed without the "openpyxl" package.');
    end
end
end