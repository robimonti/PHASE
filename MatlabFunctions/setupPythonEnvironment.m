function setupPythonEnvironment()
% setupPythonEnvironment: Automated Python Environment Setup for MATLAB

    % Check the current Python environment
    pyInfo = pyenv();
    
    % If the environment is already loaded and functional
    if pyInfo.Status == "Loaded"
        disp('Using the current Python environment:');
        disp(pyInfo.Executable);
    else
        % Prompt the user to provide the Python path
        fprintf('MATLAB is not currently using a valid Python environment.\n');
        pythonPath = input('Please provide the full path to your Python executable (e.g., /path/to/python): ', 's');
        
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

    % Share resolved Python executable with StaMPS on Windows via
    % %APPDATA%\PHASE\python.txt so the .bat shim can honor the same
    % interpreter MATLAB/PHASE is using.
    pythonPath = char(pyenv().Executable);
    if ispc
        appdata = getenv('APPDATA');
        if ~isempty(appdata)
            phase_dir = fullfile(appdata, 'PHASE');
            if ~exist(phase_dir, 'dir'); mkdir(phase_dir); end
            fid = fopen(fullfile(phase_dir, 'python.txt'), 'w', 'native', 'UTF-8');
            if fid ~= -1
                fprintf(fid, '%s', pythonPath);
                fclose(fid);
            end
        end
    end
end