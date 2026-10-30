function [lonlatIN_AOI_DET2D, dates_full, t_full, final_signal_orig, final_signal_out, final_std_out, displEXTR] = ...
    STmodel_DET2D_custom(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, x_grid, y_grid, figsDir, minMonths, ...
    gS_input_path, gS_output_path, gS_job_path, detectedOS, varNoise_method, varNoise_manual, num_spl_method, ...
    spline_method, num_spl_row_manual, num_spl_col_manual, num_spl_time_manual, lambda_method, lambda_manual, ...
    utmZone, xyIN_AOI, xyAOI, step_t, detrend_method, use_inclined_means, poly_degree, markerSize_DET2D, xy_EXTR, customConfig, customOutputDir)

% STmodel_DET1D Performs spatio-temporal deterministic modelling for displacement data
%
% Inputs:
%   displIN_AOI           - n x m matrix of displacement data (n points, m epochs)
%   PSidIN_AOI            - n x 1 vector of point IDs
%   t_dateIN              - m x 1 vector of absolute dates (datetime)
%   t_relIN               - m x 1 vector of relative times (days)
%   x_grid                - matrix containing x-coords at all query points
%   y_grid                - matrix containing y-coords at all query points
%   figsDir               - Directory path for saving figures
%   file_out              - Minimum period length for splines interpolation
%   gS_input_path         - path to input folder for geoSplinter
%   gS_output_path        - path to output folder for geoSplinter
%   gS_job_path           - path to job folder for geoSplinter
%   detectedOS            - detected operative system
%   varNoise_method       - a-priori noise variance method 'auto'/'manual'
%   varNoise_manual       - a-priori noise variance value for manual selection (otherwise NaN)
%   num_spl_method        - method for splines number selection 'auto'/'manual'
%   spline_method         - method for 'auto' splines selection 'variance'/'MDL'/'F-test'/'chi2_test'
%   num_spl_row_manual    - manual number of splines for rows (otherwise NaN)
%   num_spl_col_manual    - manual number of splines for columns (otherwise NaN)
%   num_spl_time_manual   - manual number of splines for time (otherwise NaN)
%   lambda_method         - method for lambda estimation 'auto'/'manual'
%   lambda_manual         - value of lambda for 'manual' case' (otherwise NaN)
%   utmZone               - UTM zone for reprojection of coordinates
%   xyIN_AOI              - coordinates of PS inside the AOI
%   xyAOI                 - coordiantes of the AOI limits
%   step_t                - time step for the estimation
%   detrend_method        - option to select the processing method: residual atmospheric errors (residualAtm) or clean data (cleanObs)
%   use_inclined_means    - flag for using flat or inclines plane for residual atmospheric removal
%   poly_degree           - maximum polynomial degree
%   markerSize_DET2D      - marker size for the GIF plot
%   xy_EXTR               - vector of query coordinates
% 
%
% Output:
%   lonlatIN_AOI_DET2D    - lon/lat coordinates of estimated centerline points
%   dates_full            - dates prediction
%   t_full                - relative time prediction
%   final_signal_orig     - modelled signal in observation coordinates
%   final_signal_out      - modelled signal in query coordinates
%   final_std_out         - uncertainty
%   displEXTR             - displacement at query coordinates


if nargin < 30
    xy_EXTR = [];
end
displEXTR = [];

disp('--- 2D Deterministic-based modelling started... ---');




% identify the master epoch (date where displacements are zero)
tolerance = 1e-6;    % small tolerance for numerical precision
master_idx = find(all(abs(displIN_AOI) < tolerance, 1));
if isempty(master_idx)
    warning('No master epoch found where displacements are nearly zero. Using first epoch.');
    master_idx = 1;
elseif length(master_idx) > 1
    % select epoch closest to mean of t_relIN
    [~, closest_idx] = min(abs(t_relIN(master_idx) - mean(t_relIN)));
    master_idx = master_idx(closest_idx);
    warning('Multiple near-zero displacement epochs found. Using epoch at t_relIN = %.2f.', t_relIN(master_idx));
end
master_t_rel = t_relIN(master_idx);




%% 1) SPLINES-BASED OUTLIER REJECTION
disp('======== Step 1 ========');
disp('Splines-based outlier rejection started...');

% minimum time interval in days
minT = minMonths * 30;  % months -> days

% enforce a maximum allowable interval to guarantee at least 3 splines
if round(t_relIN(end) / minT, 0) < 3
    safe_minT = t_relIN(end) / 3;
    warning('The specified minimum spline period (%.1f months) is too large for the time series duration (%.0f days). Automatically reducing it to %.1f months.', minMonths, t_relIN(end), safe_minT / 30);
    minT = safe_minT;
end

% define basis functions from geoSplinter
phi_3 = @(csi) (csi >= 0 & csi < 1) .* ((2 - csi).^3 - 4 * (1 - csi).^3) / 6 + ...
               (csi >= 1 & csi < 2) .* ((2 - csi).^3) / 6;
phi_4 = @(csi) (csi >= 0 & csi < 2) .* ((2 - csi).^3) / 6;

% initialize result storage variables
data_reg_spl = zeros(size(displIN_AOI'));
t_reg_spl = zeros(size(displIN_AOI));
data_raw_cleanT = zeros(size(displIN_AOI));

% define folder for geoSplinter based on current OS
switch detectedOS
    case 'windows'
        % Windows case
        gS_dir = (strcat('.', filesep, 'geoSplinter', filesep, 'windows')); 

    case 'linux'
        % Linux case
        gS_dir = (strcat('.', filesep, 'geoSplinter', filesep, 'linux')); 

    case 'macos_intel'
        % macOS (Intel) case
        gS_dir = (strcat('.', filesep, 'geoSplinter', filesep, 'macos_intel')); 

    case 'macos_apple_silicon'
        % macOS (Apple Silicon) case
        gS_dir = (strcat('.', filesep, 'geoSplinter', filesep, 'macos_apple_silicon')); 

    otherwise
        error('Unsupported operating system.');
end

% process each PS individually
for id = 1:size(displIN_AOI,1)
    
        fprintf('Currently processing PS %i\n', PSidIN_AOI(id));  % display current PS ID
        
        % set data and time variables for current PS
        data_reg = displIN_AOI';                      % full dataset, transposed
        t_reg = t_relIN';                              % time vector
        data_reg_PS = data_reg(:,id);                 % data for current PS
        
        % write input data to file for geoSplinter
        gS_filename = strcat('PS_', num2str(PSidIN_AOI(id)));
        writematrix([t_reg, data_reg_PS], fullfile(gS_input_path, [gS_filename, '.txt']), 'Delimiter', 'space');
        ext = 'cub';                                  % extension for cubic splines

        % define job file parameters
        data_dim = 1;                                 % data dimension
        type_spl = 2;                                 % spline type (1 = linear, 2 = cubic)
        file_inp = strcat(gS_filename, '.txt');       % input filename
        file_out = strcat(gS_filename, '_', ext);     % output filename
        num_obs = size(data_reg_PS,1);                % number of observations
        t_in = t_reg(1);                              % start time
        t_fin = round(t_reg(end));                    % end time
        lambda = 0;                                   % regularization parameter
        num_sig = 10;                                 % number of significant digits
        
        % define MDL-based spline search parameters
        max_n_spl = max(3, round(t_reg(end) / minT, 0));
        min_n_spl = 3;
        num_spl_candidates = max_n_spl - min_n_spl + 1; 

        % initialize variables
        data_spline_tmp = NaN(size(displIN_AOI, 2), num_spl_candidates);
        MDL_tmp = NaN(num_spl_candidates, 1);
    
        % initialize MDL optimization
        k = 1;                                        % iteration counter for MDL search
        for i = min_n_spl:max_n_spl
    
            num_spl = i;                              % candidate number of splines

            % write the job file
            jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, num_spl, ...
                    t_in, t_fin, lambda, num_sig, gS_input_path, gS_output_path, gS_job_path);
    
            % run geoSplinter_analysis with the job file
            gS_exec = fullfile(gS_dir, 'geoSplinter_analysis');
            gS_job_file = fullfile(gS_job_path, [file_out, '.job']);
            [status, geoOutput] = runGeoSplinter(gS_exec,gS_job_file);
            if status ~= 0
                error('Error executing geoSplinter_analysis for file %s: %s', ...
                    file_out,strtrim(geoOutput));
            end

            % define the output file base name from geoSplinter run
            mat_file = fullfile(gS_output_path, strcat(file_out, '.mat.txt'));
            hdr_file = fullfile(gS_output_path, strcat(file_out, '.hdr.txt'));

            % import normal matrix and normal vector
            fid = fopen(mat_file, 'r');
            if fid == -1
                error('Cannot open file: %s', mat_file);
            end

            n_params = num_spl + 2; % total parameters (splines + 2)
            normal_matrix = zeros(n_params, n_params);
            normal_vector = zeros(n_params, 1);
            
            while ~feof(fid)
                line = fgetl(fid);
                if ~ischar(line) || isempty(strtrim(line))
                    continue;
                end
                
                % extract Normal Matrix (not "as stored")
                if contains(line, 'Normal Matrix') && ~contains(line, 'as stored') && ~contains(line, 'Diagonal')
                    r = 1;
                    while r <= n_params
                        line = fgetl(fid);
                        if ~ischar(line) || isempty(strtrim(line))
                            continue;
                        end
                        % parse decimal format
                        values = sscanf(line, repmat('%f ', 1, n_params));
                        if length(values) == n_params
                            normal_matrix(r, :) = values';
                            r = r + 1;
                        else
                            error('Unexpected format in Normal Matrix row %d: expected %d values, got %d', ...
                                  r, n_params, length(values));
                        end
                    end
                end
                
                % extract Normal Vector
                if contains(line, 'Normal Vector')
                    r = 1;
                    while r <= n_params
                        line = fgetl(fid);
                        if ~ischar(line) || isempty(strtrim(line))
                            continue;
                        end
                        normal_vector(r) = str2double(line);
                        r = r + 1;
                    end
                end

                if contains(line, 'Sigma Zero')
                    tokens = regexp(line, 'Sigma Zero\s*:\s*([\d.]+)', 'tokens');
                    if ~isempty(tokens)
                        sigma0_est = str2double(tokens{1}{1});
                    end
                    break;
                end
            end
            fclose(fid);
            
            % import knots
            fid = fopen(hdr_file, 'r');
            if fid == -1
                error('Cannot open file: %s', hdr_file);
            end
            
            xMin = []; xMax = []; deltaX = [];
            while ~feof(fid)
                line = fgetl(fid);
                if ~ischar(line) || isempty(strtrim(line))
                    continue;
                end
                if contains(line, 'xMin')
                    xMin = sscanf(line, 'xMin   :   %e');
                elseif contains(line, 'xMax')
                    xMax = sscanf(line, 'xMax   :   %e');
                elseif contains(line, 'deltaX')
                    deltaX = sscanf(line, 'deltaX :   %e');
                end
            end
            fclose(fid);
            
            % construct knots
            knots = xMin:deltaX:xMax;

            % apply regularization for master epoch constraint
            % build design matrix A
            A = zeros(num_obs, n_params);
            for n = 1:num_obs
                t = t_reg(n);
                i_i = find(t >= knots(1:end-1) & t < knots(2:end), 1) - 1;
                if isempty(i_i)
                    if t <= knots(1)
                        i_i = 0;
                    elseif t >= knots(end)
                        i_i = length(knots) - 2;
                    end
                end
                csi_i = (t - knots(i_i + 1)) / deltaX;
                
                alpha = zeros(4, 1);
                alpha(1) = phi_4(csi_i + 1);  % B_{i_i-1}
                alpha(2) = phi_3(csi_i);      % B_{i_i}
                alpha(3) = phi_3(1 - csi_i);  % B_{i_i+1}
                alpha(4) = phi_4(2 - csi_i);  % B_{i_i+2}
                
                for s = -1:2
                    j = i_i + s + 1;
                    if j >= 1 && j <= n_params
                        A(n, j) = alpha(s + 2);
                    end
                end
            end
            
            % set B from A
            B = A(master_idx, :);
            
            % regularization
            R = B' * B;
            lambdaR = 10;
            N_modified = normal_matrix + lambdaR * R;
            x_constrained = N_modified \ normal_vector;
            
            % compute estimated observations
            estimated_obs = A * x_constrained;
            clear normal_matrix normal_vector;
            
            % calculate MDL for current spline number
            SSR = sum((data_reg_PS - estimated_obs).^2);           % sum of squared residuals
            data_spline_tmp(:,k) = estimated_obs;                  % temporary storage for current spline solution
            MDL_tmp(k) = log(num_obs) * num_spl / num_obs + log(SSR);  % MDL calculation
    
            k = k + 1;  % update iteration counter
    
        end
    
        % identify best num_spl that minimizes MDL index
        [~, MDL_idx] = min(MDL_tmp);
        
        % store final spline results and residuals
        data_spline = data_spline_tmp(:,MDL_idx);             % best spline fit
        res_det = data_reg_PS - data_spline;                  % residuals for the best spline n°
        clear data_spline_tmp MDL_tmp;

        Q = eye(size(res_det,2),size(res_det,2));
        N = A' * (Q \ A);

        % outlier removal
        Z_lim = 1.96;
        Z_obs = abs(res_det ./ sqrt(diag(sigma0_est * (Q - A / N * A'))));
        clear Q N;
        data_reg_PS(Z_obs > Z_lim) = NaN;
        t_reg(Z_obs > Z_lim) = NaN;

        % store cleaned data and time vector for this PS
        data_raw_cleanT(id,:) = data_reg_PS';
        t_reg_spl(id,:) = t_reg';
        data_reg_spl(:,id) = data_spline;
        
        fclose('all');
        pause(0.01);
end

disp('Splines-based outlier rejection completed.');




%% 2) CUSTOM TWO-STAGE MODELLING
% The preceding filter is copied verbatim from the corresponding auto engine.
% No custom mean is evaluated before this point.
[lonlatCustom, dates_full, t_full, final_signal_orig, final_signal_out, final_std_out, displEXTR] = ...
    phase_model_custom.run(data_raw_cleanT, displIN_AOI, xyIN_AOI, ...
        t_dateIN, t_relIN, [], x_grid, y_grid, xyAOI, xy_EXTR, ...
        utmZone, customConfig, customOutputDir, data_reg_spl);
lonlatIN_AOI_DET2D = lonlatCustom;
end
