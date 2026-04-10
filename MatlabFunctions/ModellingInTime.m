function [obs_p1, obs_p2, obs_p3, obs_p4, obs_p5] = ModellingInTime(...
    detectedOS, outputDir, figsDir, dataIN_AOI, displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, xyIN_AOI, ...
    varNoise_method, num_spl_method, lambda_method, coll_proc, OptionalArgs)

% ModellingInTime: performs time series interpolation.
%
% Inputs:
%   detectedOS       - Detected OS, needed for geoSplinter.
%   outputDir        - Base output directory.
%   figsDir          - Figures directory.
%   dataIN_AOI       - Matrix with PS IDs (col 1), coordinates (cols 2:3), avg velocity (col 4) and displacement time series (cols 5:end).
%   displIN_AOI      - Columns 5:end of dataIN_AOI.
%   PSidIN_AOI       - Column 1 of dataIN_AOI.
%   t_dateIN         - Vector of absolute observation dates.
%   t_relIN          - Vector of relative times (days).
%   xyIN_AOI         - Matrix of PS coordinates [x, y].
%   varNoise_method  - Noise variance method ('coherence', 'auto', 'manual').
%   num_spl_method   - Spline number selection method ('auto', 'manual').
%   lambda_method    - Regularization parameter method ('auto', 'manual').
%   coll_proc        - Collocation approach ('filtering', 'prediction').
%   OptionalArgs     - Cell array of name-value pairs for optional parameters.
%
% Optional Inputs (via OptionalArgs as name-value pairs):
%   varNoise_manual  - Manual noise variance in mm^2 (required if varNoise_method = 'manual', default: 2).
%   utmZone          - UTM zone, e.g., '33 N' (required if varNoise_method = 'coherence').
%   coherence_dir    - Directory for coherence TIFF files (required if varNoise_method = 'coherence').
%   constellation    - SAR constellation ('Sentinel-1', 'COSMO-SkyMed') (required if varNoise_method = 'coherence').
%   num_looks        - Number of looks for coherence-based variance (default: 20, optional for varNoise_method = 'coherence').
%   num_spl_manual   - Manual number of splines (required if num_spl_method = 'manual', default: 2).
%   spline_method    - Auto spline selection method ('variance', 'MDL', 'F_test', 'chi2_test', required if num_spl_method = 'auto', default: 'variance').
%   lambda_manual    - Manual regularization parameter (required if lambda_method = 'manual', default: 1).
%   coll_step_est    - Time step for prediction mode in days (required if coll_proc = 'prediction', default: 6).
%
% Outputs:
%   obs_p1 - Cell array (n_PS x 7) with detrending and outlier rejection results:
%            - Column 1: Observation dates (vector, datetime).
%            - Column 2: Relative times (vector, double, days).
%            - Column 3: Detrended displacements (vector, double, mm).
%            - Column 4: Estimated trend (vector, double, mm).
%            - Column 5: Standard deviations of polynomial coefficients (vector, double, mm).
%            - Column 6: Polynomial coefficients [a0, a1, a2] (vector, double).
%            - Column 7: Model type ('no_model', 'constant', 'linear', 'quadratic') (string).
%   obs_p2 - Cell array (n_PS x 7) with Fourier spectrum analysis results:
%            - Column 1: Frequencies (vector, double, cycles/day).
%            - Column 2: Fundamental frequencies (vector, double, cycles/day).
%            - Column 3: Noise spectrum amplitude (scalar, double, mm).
%            - Column 4: Regularized relative times (vector, double, days).
%            - Column 5: Detrended signal minus periodic component (vector, double, mm).
%            - Column 6: Fourier coefficients (vector, double).
%            - Column 7: Periodic signal (vector, double, mm).
%   obs_p3 - Cell array (n_PS x 7) with spline interpolation results:
%            - Column 1: Number of splines (scalar, integer).
%            - Column 2: Regularization parameter (scalar, double).
%            - Column 3: Observation dates (vector, datetime).
%            - Column 4: Spline-interpolated signal (vector, double, mm).
%            - Column 5: Spline residuals (vector, double, mm).
%            - Column 6: Standard deviations of spline parameters (vector, double).
%            - Column 7: Processing method description (string).
%   obs_p4 - Cell array (n_PS x 3) with covariance modeling results:
%            - Column 1: Covariance function handle (function handle).
%            - Column 2: Covariance parameters (vector, double).
%            - Column 3: Noise variance (scalar, double, mm^2).
%   obs_p5 - Cell array (n_PS x 6) with collocation and final signal reconstruction:
%            - Column 1: Final relative times (vector, double, days).
%            - Column 2: Final observation dates (vector, datetime).
%            - Column 3: Trend signal and std [signal; std] (2 x n_time, double, mm).
%            - Column 4: Periodic signal and std [signal; std] (2 x n_time, double, mm).
%            - Column 5: Collocation signal and std [signal; std] (2 x n_time, double, mm).
%            - Column 6: Final signal and std [signal; std] (2 x n_time, double, mm).


% input parser setup
p = inputParser;

% required inputs with minimal validation
addRequired(p, 'detectedOS', @ischar);
addRequired(p, 'outputDir', @ischar);
addRequired(p, 'figsDir', @ischar);
addRequired(p, 'dataIN_AOI', @ismatrix);
addRequired(p, 'displIN_AOI', @ismatrix);
addRequired(p, 'PSidIN_AOI', @isvector);
addRequired(p, 't_dateIN', @isvector);
addRequired(p, 't_relIN', @isvector);
addRequired(p, 'xyIN_AOI', @ismatrix);
addRequired(p, 'varNoise_method', @(x) ischar(x) && ismember(x, {'coherence', 'auto', 'manual'}));
addRequired(p, 'num_spl_method', @(x) ischar(x) && ismember(x, {'auto', 'manual'}));
addRequired(p, 'lambda_method', @(x) ischar(x) && ismember(x, {'auto', 'manual'}));
addRequired(p, 'coll_proc', @(x) ischar(x) && ismember(x, {'filtering', 'prediction'}));
addOptional(p, 'OptionalArgs', {}, @iscell);

% optional inputs
addParameter(p, 'varNoise_manual', 2, @(x) isscalar(x) && x > 0);
addParameter(p, 'utmZone', '', @ischar);
addParameter(p, 'coherence_dir', '', @ischar);
addParameter(p, 'constellation', 'Sentinel-1', @(x) ischar(x) && ismember(x, {'Sentinel-1', 'COSMO-SkyMed'}));
addParameter(p, 'num_looks', 20, @(x) isscalar(x) && x > 0);
addParameter(p, 'num_spl_manual', 2, @(x) isscalar(x) && x >= 2 && mod(x, 1) == 0);
addParameter(p, 'spline_method', 'variance', @(x) ischar(x) && ismember(x, {'variance', 'MDL', 'F_test', 'chi2_test'}));
addParameter(p, 'lambda_manual', 1, @(x) isscalar(x) && x > 0);
addParameter(p, 'coll_step_est', 6, @(x) isscalar(x) && x > 0);

% parse inputs
parse(p, detectedOS, outputDir, figsDir, dataIN_AOI, displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, xyIN_AOI, ...
    varNoise_method, num_spl_method, lambda_method, coll_proc, OptionalArgs);

if size(displIN_AOI, 1) ~= size(PSidIN_AOI, 1) || size(displIN_AOI, 2) ~= length(t_dateIN)
    error('Dimension mismatch: displIN_AOI, PSidIN_AOI, and t_dateIN must have consistent sizes.');
end

% override default parameters with OptionalArgs if provided
if ~isempty(OptionalArgs)
    if mod(length(OptionalArgs), 2) ~= 0
        error('OptionalArgs must contain an even number of elements (name-value pairs).');
    end
    % validate that all parameter names in OptionalArgs are recognized
    validParams = {'varNoise_manual', 'utmZone', 'coherence_dir', 'constellation', 'num_looks', ...
                   'num_spl_manual', 'spline_method', 'lambda_manual', 'coll_step_est'};
    for i = 1:2:length(OptionalArgs)
        if ~ismember(OptionalArgs{i}, validParams)
            error('Unrecognized parameter name in OptionalArgs: %s', OptionalArgs{i});
        end
    end
    % re-parse with OptionalArgs as name-value pairs
    parse(p, detectedOS, outputDir, figsDir, dataIN_AOI, displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, xyIN_AOI, ...
        varNoise_method, num_spl_method, lambda_method, coll_proc, OptionalArgs, OptionalArgs{:});
end

% extract optional inputs
varNoise_manual = p.Results.varNoise_manual;
utmZone = p.Results.utmZone;
coherence_dir = p.Results.coherence_dir;
constellation = p.Results.constellation;
num_looks = p.Results.num_looks;
num_spl_manual = p.Results.num_spl_manual;
spline_method = p.Results.spline_method;
lambda_manual = p.Results.lambda_manual;
coll_step_est = p.Results.coll_step_est;

% conditional input validation
% - varNoise_method
switch varNoise_method
    case 'coherence'
        if isempty(utmZone)
            error('utmZone must be provided when varNoise_method is ''coherence''.');
        end
        if isempty(coherence_dir)
            error('coherence_dir must be provided when varNoise_method is ''coherence''.');
        end
        if isempty(constellation)
            error('constellation must be provided when varNoise_method is ''coherence''.');
        end
        if ~exist(coherence_dir, 'dir')
            error('coherence_dir ''%s'' does not exist.', coherence_dir);
        end
    case 'manual'
        if isnan(varNoise_manual) || varNoise_manual <= 0
            error('varNoise_manual must be a positive scalar when varNoise_method is ''manual''.');
        end
end
% - num_spl_method
switch num_spl_method
    case 'auto'
        if ~ismember(spline_method, {'variance', 'MDL', 'F_test', 'chi2_test'})
            error('spline_method must be one of ''variance'', ''MDL'', ''F_test'', or ''chi2_test'' when num_spl_method is ''auto''.');
        end
    case 'manual'
        if isnan(num_spl_manual) || num_spl_manual < 2 || mod(num_spl_manual, 1) ~= 0
            error('num_spl_manual must be an integer >= 2 when num_spl_method is ''manual''.');
        end
end
% - lambda_method
switch lambda_method
    case 'manual'
        if isnan(lambda_manual) || lambda_manual <= 0
            error('lambda_manual must be a positive scalar when lambda_method is ''manual''.');
        end
end
% - coll_proc
switch coll_proc
    case 'prediction'
    if isnan(coll_step_est) || coll_step_est <= 0
        error('coll_step_est must be a positive scalar when coll_proc is ''prediction''.');
    end
end

% set lambda_sar based on constellation
lambda_s1 = 0.05546576;   % Sentinel-1 wavelength [m]
lambda_csk = 0.031228;    % COSMO-SkyMed wavelength [m]
switch constellation
    case 'Sentinel-1'
        lambda_sar = lambda_s1;
    case 'COSMO-SkyMed'
        lambda_sar = lambda_csk;
    otherwise
        error('Invalid constellation: %s. Must be ''Sentinel-1'' or ''COSMO-SkyMed''.', constellation);
end

% compute coll_t_est for prediction mode
if strcmp(coll_proc, 'prediction')
    if isempty(t_relIN)
        error('t_relIN is empty, cannot define coll_t_est.');
    end
    t_relIN_sorted = sort(t_relIN);
    coll_t_est = t_relIN_sorted(1):coll_step_est:t_relIN_sorted(end);
    if isempty(coll_t_est)
        error('coll_t_est is empty for PS %d. Check coll_step_est and t_relIN.', PSidIN_AOI(i));
    end
end


% check for required functions in MatlabFunctions folder
requiredFunctions = {'read_coherence_tifs', 'compute_varNoise_from_coherence1D', ...
                    'jobFile_analysis', 'geoSplinter_noFig', 'jobFile_synthesis', ...
                    'compute_B_spline_row', 'f1DEmpCovEst', 'nearestSPD'};
for func = requiredFunctions
    if ~exist(fullfile('MatlabFunctions', func{1}), 'file')
        error('Required function %s not found in MatlabFunctions folder.', func{1});
    end
end

% ----------- end of input variables preparation -----------


disp('--- PS time series modelling started... ---');


% create processing folders for geoSplinter
gS_base_path = fullfile(outputDir, 'geoSplinter');
folders = {'data_input', 'data_output', 'job', 'data_synthesis'};
for folder = folders
    dir_path = fullfile(gS_base_path, folder{1});
    if ~exist(dir_path, 'dir')
        [success, msg] = mkdir(dir_path);
        if ~success
            error('Failed to create directory %s: %s', dir_path, msg);
        end
    end
end
gS_input_path = fullfile(gS_base_path, 'data_input');
gS_output_path = fullfile(gS_base_path, 'data_output');
gS_job_path = fullfile(gS_base_path, 'job');
gS_synth_path = fullfile(gS_base_path, 'data_synthesis');


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
master_date = t_dateIN(master_idx);
master_t_rel = t_relIN(master_idx);

% initialize variables
obs_p1 = cell(size(displIN_AOI, 1), 6);
obs_p2 = cell(size(displIN_AOI, 1), 7);
obs_p3 = cell(size(displIN_AOI, 1), 7);
obs_p4 = cell(size(displIN_AOI, 1), 3);
obs_p5 = cell(size(displIN_AOI, 1), 6);

% maximum number of iterations
max_iterations = 10;

for i = 1:size(dataIN_AOI, 1)

    % counter for iterations
    ii = 1;

    fprintf('======== PS %i ========\n', PSidIN_AOI(i));
    disp('Processing started ...')

    obs = displIN_AOI(i,:);
    t0_date = t_dateIN(1);

    % -- 1) Detrending and outlier rejection
    % normalize time values to improve numerical stability
    t_tmp = t_relIN;
    t_mean = mean(t_relIN);
    t_std = std(t_relIN);
    t_norm = (t_relIN - t_mean) / t_std;
    
    % initial design matrix with normalized time values
    A = [ones(length(t_norm), 1), t_norm', t_norm'.^2];
    obs_col = obs(:);
    alpha_sig = 0.05;
    dof = size(A,1) - size(A,2);

    % t-test: 0 = constant, 1 = linear, 2 = quadratic
    terms = [0, 1, 2];

    % create weight matrix (inverse of Q matrix)
    W = eye(length(obs_col));
    master_idx_tmp = find(abs(t_tmp - master_t_rel) < eps, 1);
    % assign very low variance (high weight) to master epoch
    if ~isempty(master_idx_tmp)
        master_variance = 1e-1;
        W(master_idx_tmp, master_idx_tmp) = 1 / master_variance;
    else
        warning('Master epoch not found in time vector for PS %d.', PSidIN_AOI(i));
    end
    
    % iterative model reduction
    while true
        % LS estimation
        if cond(A) > 1e10
            x_est = pinv(A' * W * A) * (A' * W * obs_col);
        else
            x_est = (A' * W * A) \ (A' * W * obs_col);
        end
        y_est = A * x_est;
        v_est = obs_col - y_est;
        s02_est = (v_est' * W * v_est) / dof;

        % covariance matrix
        ATA = A' * W * A;
        Cxx_est = s02_est * (ATA \ eye(size(ATA)));
    
        % compute t-values
        t_x_est = x_est ./ sqrt(diag(Cxx_est));
        t_lim = tinv(1 - alpha_sig/2, dof);
    
        % find the least significant parameter
        [min_t_val, min_idx] = min(abs(t_x_est));
    
        if min_t_val >= t_lim
            % stop if all parameters are significant
            break;
        end
    
        % remove the least significant parameter
        A(:, min_idx) = [];
        terms(min_idx) = [];
        dof = size(A,1) - size(A,2);

        % check if all parameters are removed
        if isempty(A) || size(A,2) == 0
            % no terms left, exit loop
            break;
        end

    end

    % OUTLIER REJECTION
    if isempty(terms)

        % if all terms were removed set coefficients to zero
        a0 = 0; a1 = 0; a2 = 0;
        y_est_final = zeros(size(t_tmp));
        std_poly_final = NaN(3, 1);
        model_type = 'no_model';

    else

        % outlier rejection
        max_outliers = floor(0.05 * length(obs_col));
        outliers_removed = 0;
        obs_out = NaN(max_outliers, 1);
        t_out = NaN(max_outliers, 1);
    
        A_tmp = A;
        obs_col_tmp = obs_col;
        t_tmp_tmp = t_tmp;
        W_tmp = W;
    
        while outliers_removed < max_outliers
            % weighted LS estimation with current data
            if cond(A_tmp) > 1e10
                x_est = pinv(A_tmp' * W_tmp * A_tmp) * (A_tmp' * W_tmp * obs_col_tmp);
            else
                x_est = (A_tmp' * W_tmp * A_tmp) \ (A_tmp' * W_tmp * obs_col_tmp);
            end
            y_est = A_tmp * x_est;
            v_est = obs_col_tmp - y_est;
    
            % calculate MAD for residuals
            median_res = median(v_est);
            mad_var = median(abs(v_est - median_res)) * 1.4826;
            if mad_var < eps, mad_var = eps; end
    
            % standardized residuals using MAD
            robust_std_res = (v_est - median_res) / mad_var;
    
            % find the largest absolute residual
            [max_abs_res, max_idx] = max(abs(robust_std_res));
            outlier_threshold = 3;
    
            if max_abs_res <= outlier_threshold
                % stop if no more outliers exceed threshold
                break;
            end

            % prevent master epoch from being flagged
            if max_idx == master_idx_tmp
                robust_std_res(max_idx) = 0; 
                [~, max_idx] = max(abs(robust_std_res));
            end
    
            % remove the outlier
            obs_out(outliers_removed + 1) = obs_col_tmp(max_idx);
            t_out(outliers_removed + 1) = t_tmp_tmp(max_idx);
            A_tmp(max_idx, :) = [];
            obs_col_tmp(max_idx) = [];
            t_tmp_tmp(max_idx) = [];
            W_tmp(max_idx, :) = [];
            W_tmp(:, max_idx) = [];
            outliers_removed = outliers_removed + 1;
        end
    
        % final weighted LS estimation after outlier rejection
        if cond(A_tmp) > 1e10
            x_est_final = pinv(A_tmp' * W_tmp * A_tmp) * (A_tmp' * W_tmp * obs_col_tmp);
        else
            x_est_final = (A_tmp' * W_tmp * A_tmp) \ (A_tmp' * W_tmp * obs_col_tmp);
        end
    
        % get the error estimate
        y_est = A_tmp * x_est_final;
        v_est = obs_col_tmp - y_est;
        s02_est = (v_est' * W_tmp * v_est) / (size(A_tmp,1) - size(A_tmp,2));
        ATA = A_tmp' * W_tmp * A_tmp;
        Cxx_est = s02_est * (ATA \ eye(size(ATA)));
        std_poly = sqrt(diag(Cxx_est));
    
        % initialize coefficients
        a0 = 0; a1 = 0; a2 = 0;
        std_poly_final = NaN(3, 1);
    
        % convert back to original scale
        % map coefficients and std based on remaining terms
        for pp = 1:length(terms)
            switch terms(pp)
                case 0 % constant term
                    a0 = x_est_final(pp);
                    std_poly_final(1) = std_poly(pp);
                case 1 % linear term
                    a1 = x_est_final(pp) / t_std;
                    a0 = a0 - x_est_final(pp) * t_mean / t_std;
                    std_poly_final(2) = std_poly(pp) / t_std;
                case 2 % quadratic term
                    a2 = x_est_final(pp) / (t_std^2);
                    a1 = a1 - 2 * x_est_final(pp) * t_mean / (t_std^2);
                    a0 = a0 + x_est_final(pp) * (t_mean^2) / (t_std^2);
                    std_poly_final(3) = std_poly(pp) / (t_std^2);
            end
        end
        
        % final estimated trend in original scale (for all original points)
        y_est_final = a0 + a1 * t_tmp + a2 * t_tmp.^2;
        
        max_degree = max(terms);
        switch max_degree
            case 0
                model_type = 'constant';
            case 1
                model_type = 'linear';
            case 2
                model_type = 'quadratic';
        end

    end
    
    % detrended original data
    obs_detrended = obs_col - y_est_final(:);
    t_date_out = t_tmp + t0_date;

    % save outlier-free observations and time vectors
    obs_p1{i, 1} = t_date_out;
    obs_p1{i, 2} = t_tmp;
    obs_p1{i, 3} = obs_detrended';
    obs_p1{i, 4} = y_est_final;
    obs_p1{i, 5} = std_poly_final;
    obs_p1{i, 6} = [a0; a1; a2];
    obs_p1{i, 7} = model_type;


    % initialize the variance of the noise
    switch varNoise_method
        case 'manual'
            if ~isnumeric(varNoise_manual) || ~isscalar(varNoise_manual) || varNoise_manual <= 0
                warning('Invalid varNoise_manual (%s) for PS %d. Using default: 2 mm^2.', ...
                    num2str(varNoise_manual), PSidIN_AOI(i));
                varNoise = 2;
            else
                varNoise = varNoise_manual;
            end
        case 'auto'
            varNoise = varNoise_manual; % initial fallback
        case 'coherence'
            if i == 1
                coh_values = read_coherence_tifs(xyIN_AOI, t_relIN - master_t_rel, ...
                    coherence_dir, master_date, utmZone);
                varNoise_all = compute_varNoise_from_coherence1D(coh_values, ...
                    lambda_sar, num_looks, varNoise_manual);
            end
            varNoise = varNoise_all(i);
            if isnan(varNoise) || varNoise <= 0
                warning('Invalid coherence-based varNoise for PS %d. Using manual: %.2f mm^2.', ...
                    PSidIN_AOI(i), varNoise_manual);
                varNoise = varNoise_manual;
            end
    end

    varNoise_initial = varNoise; % store initial value for fallback
    varNoise_cov_prec = []; % initialize as empty to indicate no previous value


    while ii <= max_iterations

        % -- 2) Fourier spectrum for outliers
        % compute the mode for uniform sampling
        mode_obs = mode(diff(obs_p1{i, 2}));
    
        % compute the regularized time vector
        t_reg = min(obs_p1{i, 2}) : mode_obs : max(max(obs_p1{i, 2}), mode_obs * ceil(max(obs_p1{i, 2})/mode_obs));
    
        % re-sample the signal at the mode with interp1
        obs_reg = interp1(obs_p1{i, 2}, obs_p1{i, 3}, t_reg, "linear");
    
        % compute the PSD and amplitude spectrum
        [obs_psd, freq, ~, ~] = sigs2epsd(obs_reg, t_reg', '--noPadding');
        obs_spectr = sqrt(obs_psd);
    
        % moving median on the spectrum
        win_size = min(0.1 * length(obs_p1{i, 2}), 50);
        win_size = max(win_size, 3); % ensure odd window size for movmedian
        mov_spectr = movmedian(obs_spectr, win_size, 'omitnan');        
    
        % residuals between amplitude spectrum and movmedian
        diff_spectr = obs_spectr - mov_spectr;

        % check for clear periodicity
        max_amplitude = max(obs_spectr);
        median_background = median(obs_spectr);
        mad_background = 1.4826 * mad(obs_spectr, 1);
        periodicity_ratio = max_amplitude / max(median_background, mad_background);
        periodicity_threshold = 10;

        % initialize master_idx_reg
        master_idx_reg = find(abs(t_reg - master_t_rel) < eps, 1);
        if isempty(master_idx_reg)
            warning('Master epoch not found in regularized time vector for PS %d.', PSidIN_AOI(i));
        end

        if periodicity_ratio > periodicity_threshold
            % CLEAR PERIODICITY EXISTS
            % determination of fundamental frequencies
            % identify outliers (potential peaks) using GESD
            outlier_idx_gesd = isoutlier(abs(diff_spectr), 'gesd');
        
            % define a threshold for identifying dominant harmonics using "outliers"
            outlier_mags = diff_spectr(outlier_idx_gesd);
            outlier_mags_full = diff_spectr;
            outlier_mags_full(~outlier_idx_gesd) = 0;
            idx_cen = ceil(length(freq) / 2);
            outlier_mags_half = outlier_mags_full(idx_cen+1:end);

            if isempty(outlier_mags)
                % fallback if no outliers detected
                fund_threshold = max(max(diff_spectr) * 0.1, eps);
                [~, fund_idx] = findpeaks(outlier_mags_half, 'MinPeakHeight', fund_threshold);
            else
                % findpeaks to identify fundamental frequencies
                nonzero_outlier_mags_half = outlier_mags_half(outlier_mags_half ~= 0);
                if isempty(nonzero_outlier_mags_half)
                    noiseLevel = eps;
                else
                    noiseLevel = 1.4826 * mad(nonzero_outlier_mags_half, 1);
                end
                minProminence = max(0.80 * max(outlier_mags_half), 2 * noiseLevel); 
    
                % calculate the minimum peak distance (5% of the half-spectrum length)
                minDistance = round(length(outlier_mags_half) / 20); 
            
                % find peaks in the "outliers" spectrum
                [~, fund_idx] = findpeaks(outlier_mags_half, 'MinPeakProminence', minProminence, 'MinPeakDistance', minDistance);
            end
        
            % identify fundamental frequencies
            fundamental_freq = freq(idx_cen + fund_idx);
        else
            % NO CLEAR PERIODICITY DETECTED
            % no fundamental frequencies are set
            fundamental_freq = [];
        end

        % reconstruct periodic signal using cosine series
        if ~isempty(fundamental_freq)
            % construct design matrix for cosine series
            n_freq = length(fundamental_freq);
            A_periodic = zeros(length(t_reg), 2 * n_freq + 1);                     % cosines, sines, and constant term
            A_periodic(:, 1) = 1;                                                  % constant term
            for k = 1:n_freq
                A_periodic(:, 2*k) = cos(2 * pi * fundamental_freq(k) * t_reg');   % cosine terms
                A_periodic(:, 2*k+1) = sin(2 * pi * fundamental_freq(k) * t_reg'); % sine terms
            end
        
            % observation vector (detrended signal)
            y_periodic = obs_reg(:);
        
            % weight matrix (high weight for master epoch)
            W_periodic = eye(length(t_reg));
            if ~isempty(master_idx_reg)
                master_variance = 1e-1; 
                W_periodic(master_idx_reg, master_idx_reg) = 1 / master_variance;
            else
                warning('Master epoch not found in regularized time vector for PS %d.', PSidIN_AOI(i));
            end
        
            % weighted LS estimation
            if cond(A_periodic) > 1e10
                coeffs = pinv(A_periodic' * W_periodic * A_periodic) * (A_periodic' * W_periodic * y_periodic);
            else
                coeffs = (A_periodic' * W_periodic * A_periodic) \ (A_periodic' * W_periodic * y_periodic);
            end
        
            % reconstruct periodic signal
            mod_time_signal = (A_periodic * coeffs)';
        else
            % no periodicity: set periodic signals to zero
            mod_time_signal = zeros(size(t_reg));
            fundamental_freq = [];
            coeffs = [];
        end

        % compute noise_spectr (and varNoise for 'auto' method) in first iteration
        if ii == 1
            if ~isempty(fundamental_freq)
                % periodic signal detected
                residuals = y_periodic' - mod_time_signal;
                noise_spectr = 1.4826 * mad(residuals, 1);
            else
                % no periodicity: use MAD of detrended observations
                noise_spectr = 1.4826 * mad(obs_reg, 1);
            end
            % handle invalid noise_spectr
            if noise_spectr <= 0 || isnan(noise_spectr)
                warning('Invalid noise_spectr (%.2f) from Fourier analysis for PS %d. Using default: %.2f mm.', ...
                    noise_spectr, PSidIN_AOI(i), sqrt(varNoise_manual));
                noise_spectr = sqrt(varNoise_manual);
            end
            % update varNoise only for 'auto' method
            switch varNoise_method
                case 'auto'
                varNoise = noise_spectr^2;
                if varNoise <= 0 || isnan(varNoise)
                    warning('Invalid varNoise (%.2f) from Fourier analysis for PS %d. Using default: %.2f mm^2.', ...
                        varNoise, PSidIN_AOI(i), varNoise_manual);
                    varNoise = varNoise_manual;
                end
            end
        else
            % use updated varNoise from covariance step
            noise_spectr = sqrt(varNoise);
        end
        
        % sample reconstructed signals at observation times
        tobs_idx = ismember(t_reg, obs_p1{i, 2});
        t_rel_orig = t_reg(tobs_idx);
        mod_time_signal_orig_fund = mod_time_signal(tobs_idx);
        mod_time_signal_orig = obs_p1{i, 3} - mod_time_signal_orig_fund;

        % save outlier free spectrum and time series
        obs_p2{i, 1} = freq';
        obs_p2{i, 2} = fundamental_freq';
        obs_p2{i, 3} = noise_spectr;
        obs_p2{i, 4} = t_rel_orig;
        obs_p2{i, 5} = mod_time_signal_orig;
        obs_p2{i, 6} = coeffs;
        obs_p2{i, 7} = mod_time_signal_orig_fund;
    
    
        % -- 3) Splines interpolation with geoSplinter
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
    
        % write the input file
        gS_filename = strcat('PS_', num2str(PSidIN_AOI(i)));
        writematrix([obs_p2{i, 4}', obs_p2{i, 5}'], fullfile(gS_input_path, [gS_filename, '.txt']), 'Delimiter', 'space')
        
        % define processing variables
        ext = 'cub';                                    % ext variable based on splines type
    
        % job file creation
        data_dim = 1;                                   % Dimension of dataset (1D)
        type_spl = 2;                                   % Type of splines (cubic)
        file_inp = strcat(gS_filename, '.txt');         % Input filename (with extension)
        file_out = strcat(gS_filename, '_', ext);       % Output filename (without extension)
        num_obs  = length(obs_p2{i, 5});                % Number of observations
        t_in     = obs_p1{i, 2}(1);                     % First abscissa (= time)
        t_fin    = obs_p1{i, 2}(end);                   % Last abscissa (= time)
        num_sig  = 10;                                  % Number of significant digits
    
        % define min and max number of splines to be tested
        min_n_spl = 2;
    
        T_min = 5 * median(diff(obs_p2{i, 4}));
        max_n_spl = min(round((max(obs_p2{i, 4}) - min(obs_p2{i, 4})) / T_min), round(num_obs * 0.1));
        max_n_spl = max(max_n_spl, 2);

        % validate manual inputs
        switch num_spl_method
            case 'manual'
            if ~isnumeric(num_spl_manual) || ~isscalar(num_spl_manual) || num_spl_manual < min_n_spl || mod(num_spl_manual, 1) ~= 0
                warning('Invalid num_spl_manual (%s) for PS %d. Reverting to auto.', num2str(num_spl_manual), PSidIN_AOI(i));
                num_spl_method = 'auto';
            end
        end
        switch lambda_method
            case 'manual'
            if ~isnumeric(lambda_manual) || ~isscalar(lambda_manual) || lambda_manual <= 0
                warning('Invalid lambda_manual (%s) for PS %d. Reverting to auto.', num2str(lambda_manual), PSidIN_AOI(i));
                lambda_method = 'auto';
            end
        end

        % define support functions for geoSplinter
        phi_3 = @(csi) (csi >= 0 & csi < 1) .* ((2 - csi).^3 - 4 * (1 - csi).^3) / 6 + ...
               (csi >= 1 & csi < 2) .* ((2 - csi).^3) / 6;
        phi_4 = @(csi) (csi >= 0 & csi < 2) .* ((2 - csi).^3) / 6;
    
        % initialize variables
        switch num_spl_method
            case 'manual'
                % only one run for manual num_spl
                num_spl_candidates = 1; 
            case 'auto'
                % multiple runs for auto
                num_spl_candidates = max_n_spl - min_n_spl + 1; 
        end
        PS_data_spline_tmp = NaN(size(obs_p2{i,5}, 2), num_spl_candidates);
        s02_spl_tmp = NaN(1, num_spl_candidates);
        var_resSpl = NaN(num_spl_candidates, 1);
        MDL_tmp = NaN(num_spl_candidates, 1);
        F_pval_tmp = NaN(num_spl_candidates, 1);
        chi2_pval_tmp = NaN(num_spl_candidates, 1);
        lambda_it = NaN(num_spl_candidates, 1);
        diag_Cxx = cell(num_spl_candidates, 1);
        diag_Cxx_full = cell(num_spl_candidates, 1);
        knots_tmp = cell(num_spl_candidates, 1);

        % determine splines number to test
        switch num_spl_method
            case 'manual'
                splines_to_test = num_spl_manual;
            case 'auto'
                splines_to_test = min_n_spl:max_n_spl;
        end
    
        % the number of splines is iteratively chosen
        kk = 1;
        for j = splines_to_test
    
            % incremental definition of num_spl
            num_spl = j;

            % number of parameters for cubic splines
            n_params = num_spl + 2;
    
            % regularization parameter (λ)
            switch lambda_method
                case 'auto'
                    tauGrid = (obs_p2{i, 4}(end) - obs_p2{i, 4}(1)) / (num_spl - 1);
                    lambda = lambdaSplines1D([obs_p2{i, 4}', obs_p2{i, 5}'], varNoise, tauGrid, type_spl);
                case 'manual'
                    lambda = lambda_manual;
            end
    
            % write the job file
            jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, num_spl, ...
                    t_in, t_fin, lambda, num_sig, gS_input_path, gS_output_path, gS_job_path);
    
            % run geoSplinter_analysis with the job file
            job_execution = sprintf('%s < %s', fullfile(gS_dir, 'geoSplinter_analysis'), fullfile('.', gS_job_path, [file_out, '.job']));
            status = system(job_execution);
            if status ~= 0
                error('Error executing geoSplinter_analysis for file: %s', file_out);
            end
            
            % import results
            [gS_out_data_tmp, ~, ~, ~, ~] = ...
                geoSplinter_noFig(fullfile('.', gS_output_path, file_out), ext);
    
            % import normal matrix and normal vector
            mat_file = fullfile('.', gS_output_path, strcat(file_out, '.mat.txt'));
            fid = fopen(mat_file, 'r');
            if fid == -1
                error('Cannot open file: %s', mat_file);
            end
            try
                normal_matrix = zeros(n_params, n_params);
                normal_vector = zeros(n_params, 1);
                while ~feof(fid)
                    line = fgetl(fid);
                    if ~ischar(line) || isempty(strtrim(line))
                        continue;
                    end
                    if contains(line, 'Normal Matrix') && ~contains(line, 'as stored') && ~contains(line, 'Diagonal')
                        r = 1;
                        while r <= n_params
                            line = fgetl(fid);
                            if ~ischar(line) || isempty(strtrim(line))
                                continue;
                            end
                            values = sscanf(line, repmat('%f ', 1, n_params));
                            if length(values) == n_params
                                normal_matrix(r, :) = values';
                                r = r + 1;
                            end
                        end
                    end
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
                end
            catch err
                fclose(fid);
                error('Error reading file %s: %s', mat_file, err.message);
            end
            fclose(fid);
        
            % import knots
            hdr_file = fullfile('.', gS_output_path, strcat(file_out, '.hdr.txt'));
            fid = fopen(hdr_file, 'r');
            if fid == -1
                error('Cannot open file: %s', hdr_file);
            end
            try
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
            catch err
                fclose(fid);
                error('Error reading file %s: %s', hdr_file, err.message);
            end
            fclose(fid);
            knots = xMin:deltaX:xMax;
            if isempty(knots) || length(knots) < 2
                error('Invalid knot vector for PS %d: knots must have at least 2 elements.', PSidIN_AOI(i));
            end
        
            % apply regularization for master epoch constraint
            master_idx = find(abs(obs_p2{i, 4} - master_t_rel) < eps, 1);
            if ~isempty(master_idx)
                B = compute_B_spline_row(master_t_rel, knots, deltaX, n_params);
                R = B' * B;
                lambdaR = 1e1; % regularization strength
                N_modified = normal_matrix + lambdaR * R;
                x_constrained = N_modified \ normal_vector;
                % compute estimated observations
                A = zeros(num_obs, n_params);
                for n = 1:num_obs
                    t = obs_p2{i, 4}(n);
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
                    alpha(1) = phi_4(csi_i + 1);
                    alpha(2) = phi_3(csi_i);
                    alpha(3) = phi_3(1 - csi_i);
                    alpha(4) = phi_4(2 - csi_i);
                    for s = -1:2
                        u = i_i + s + 1;
                        if u >= 1 && u <= n_params
                            A(n, u) = alpha(s + 2);
                        end
                    end
                end
                estimated_obs = A * x_constrained;
            else
                warning('Master epoch not found in time vector for PS %d.', PSidIN_AOI(i));
                estimated_obs = gS_out_data_tmp(:, 3); % fallback to original geoSplinter output
            end

            % compute residual variance and SSR
            residuals = obs_p2{i, 5} - estimated_obs';
            SSR = sum(residuals.^2);
            s02_est = SSR / (num_obs - n_params);
            gS_out_mat = readmatrix(fullfile('.', gS_output_path, strcat(file_out, '.mat.txt')));
            gS_out_mat_col = mean(gS_out_mat, 2, 'omitnan');
            idx_lnan = find(isnan(gS_out_mat_col), 1, 'last');
            diag_Cxx{kk} = sqrt(gS_out_mat_col(idx_lnan+2:end-2));
            diag_Cxx_full{kk} = sqrt(gS_out_mat_col(idx_lnan+1:end-1));
            knots_tmp{kk} = knots;

            switch num_spl_method
                case 'auto'
                    % compute MDL index
                    MDL_tmp(kk) = log(num_obs) * num_spl / num_obs + log(max(SSR, eps));
                
                    % compute chi2 test p-value
                    if num_obs > n_params && varNoise > 0
                        chi2_stat = SSR / varNoise;
                        dof = num_obs - n_params;
                        chi2_pval_tmp(kk) = 1 - chi2cdf(chi2_stat, dof);
                    else
                        chi2_pval_tmp(kk) = NaN; 
                    end
            end

            % store splines interpolations
            PS_data_spline_tmp(:,kk) = estimated_obs;
            s02_spl_tmp(kk) = s02_est;
            var_resSpl(kk) = s02_est;
            lambda_it(kk) = lambda;
        
            % update the counter
            kk = kk + 1;
    
        end

        % compute F-test p-values (compare num_spl vs. num_spl-1)
        switch num_spl_method
            case 'auto'
                for kk = 2:num_spl_candidates
                    num_spl = min_n_spl + kk - 1;
                    n_params = num_spl + 2;
                    dof_complex = num_obs - n_params;
                    SSR_complex = s02_spl_tmp(kk) * dof_complex;
                    
                    num_spl_simple = num_spl - 1;
                    n_params_simple = num_spl_simple + 2;
                    dof_simple = num_obs - n_params_simple;
                    SSR_simple = s02_spl_tmp(kk-1) * dof_simple;
                    
                    if SSR_simple > SSR_complex && dof_complex > 0
                        F_stat = ((SSR_simple - SSR_complex) / (dof_simple - dof_complex)) / s02_spl_tmp(kk);
                        F_pval_tmp(kk) = 1 - fcdf(F_stat, dof_simple - dof_complex, dof_complex);
                    else
                        F_pval_tmp(kk) = NaN;
                    end
                end
                F_pval_tmp(1) = 0;
        end
        
        % select number of splines
        switch num_spl_method
            case 'manual'
                num_spl = num_spl_manual;
                var_idx = 1;
            case 'auto'
                % select number of splines based on spline_method for 'auto'
                switch spline_method
                    case 'variance'
                        % variance-based method
                        var_thrs = prctile(var_resSpl, 80);
                        var_spl_diff = var_resSpl - var_thrs;
                        var_spl_diff(var_spl_diff < 0) = nan;
                        [~, var_idx] = min(var_spl_diff, [], 'omitnan');
                        if ~isempty(var_idx)
                            num_spl = min_n_spl + var_idx - 1;
                        else
                            num_spl = min_n_spl;
                            var_idx = 1;
                        end
                    case 'MDL'
                        % MDL-based method
                        [~, MDL_idx] = min(MDL_tmp);
                        num_spl = min_n_spl + MDL_idx - 1;
                        var_idx = MDL_idx;
                    case 'F_test'
                        % F-test method
                        if all(isnan(F_pval_tmp))
                            warning('F-test p-values are all NaN for PS %d, defaulting to min_n_spl.', PSidIN_AOI(i));
                            num_spl = min_n_spl;
                            var_idx = 1;
                        else
                            % select smallest num_spl where p-value > alpha (no significant improvement)
                            idx = find(F_pval_tmp > alpha_sig, 1, 'first');
                            if isempty(idx)
                                % if all p-values are significant, choose max_n_spl
                                num_spl = max_n_spl;
                                var_idx = num_spl_candidates;
                            else
                                % select num_spl-1 (simpler model)
                                num_spl = min_n_spl + idx - 2;
                                var_idx = idx - 1;
                                if num_spl < min_n_spl
                                    num_spl = min_n_spl;
                                    var_idx = 1;
                                end
                            end
                        end
                    case 'chi2_test'
                        % chi-squared test method
                        if all(isnan(chi2_pval_tmp))
                            warning('Chi2 test p-values are all NaN for PS %d, defaulting to min_n_spl.', PSidIN_AOI(i));
                            num_spl = min_n_spl;
                            var_idx = 1;
                        else
                            % select smallest num_spl where p-value is in acceptable range
                            chi2_pval_range = [alpha_sig/2, 1-alpha_sig/2];
                            idx = find(chi2_pval_tmp >= chi2_pval_range(1) & chi2_pval_tmp <= chi2_pval_range(2), 1, 'first');
                            if isempty(idx)
                                % if no p-value is in range, choose min_n_spl to avoid overfitting
                                num_spl = min_n_spl;
                                var_idx = 1;
                            else
                                num_spl = min_n_spl + idx - 1;
                                var_idx = idx;
                            end
                        end
                    otherwise
                        error('Invalid spline_method: %s. Choose ''variance'', ''MDL'', ''F_test'', ''chi2_test''.', spline_method);
                end
        end

        % get the result from var_idx
        PS_data_spline = PS_data_spline_tmp(:, var_idx);
        res_spl = obs_p2{i, 5} - PS_data_spline';
        lambda_spl = lambda_it(var_idx);
        std_splParms = diag_Cxx{var_idx};
        std_splParms_full = diag_Cxx_full{var_idx};
        knots_csn = knots_tmp{var_idx};
        
        % populate output matrices
        obs_p3{i, 1} = num_spl;
        obs_p3{i, 2} = lambda_spl;
        obs_p3{i, 3} = t0_date + obs_p2{i, 4};
        obs_p3{i, 4} = PS_data_spline';
        obs_p3{i, 5} = res_spl;
        obs_p3{i, 6} = std_splParms;
        switch num_spl_method
            case 'auto'
                obs_p3{i, 7} = sprintf('num_spl:auto(%s),lambda:%s,varNoise:%s', spline_method, lambda_method, varNoise_method);
            case 'manual'
                obs_p3{i, 7} = sprintf('num_spl:manual,lambda:%s,varNoise:%s', lambda_method, varNoise_method);
        end
    
        clear PS_data_spline
        clear res_spl
    
        
        % -- 4) Covariance modelling
        % define the possible empirical models
        % gaussian
        mCovF1 = @(c, tau) c(1) .* exp(-c(2) .* tau.^2);
        % gaussian with bell
        mCovF2 = @(c, tau) c(1) .* exp(-c(2) .* tau.^2) .* (1 - c(3) .* tau .^2);
        % gaussian with cosine
        mCovF3 = @(c, tau) c(1) .* exp(-c(2) .* tau.^2) .* cos(c(3) .* tau .* pi/2);
        
        % compute the empirical covariance function
        dtau = 2 * mode_obs;
        [~, ~, tauGrid, eCovF, Cecf, ~] = f1DEmpCovEst(obs_p3{i, 5}', obs_p2{i, 4}', dtau, 0);
        
        % create temporary copies for polynomial fitting
        tauGrid_poly = tauGrid;
        eCovF_poly = eCovF;
        Cecf_poly = Cecf;
        
        % check if second value is negative and add constraint in temporary copies
        if eCovF(2) < 0
            x = [tauGrid_poly(2), tauGrid_poly(3)];
            y = [eCovF_poly(2), eCovF_poly(3)];
            tau0_interp = interp1(x, y, 0, 'linear', 'extrap');
            tau0_interp = max(tau0_interp, 0);
            tauGrid_poly = [tauGrid_poly(1); 0; tauGrid_poly(2:end)];
            eCovF_poly = [eCovF_poly(1); tau0_interp; eCovF_poly(2:end)];
            Cecf_poly = [Cecf_poly(1); Cecf_poly(1); Cecf_poly(2:end)];
        end
        
        % fit each empirical covariance model
        % find the tau coordinates for which the covariance function becomes 0
        idxZero = find(eCovF < 0, 1, 'first');
        
        % extract covariance points to be used in interpolation
        idx_in = 2;
        idx_end_min = 15;    % ensure minimum number of points
        idx_end = max(idx_end_min, round(size(Cecf,1)/2));
        
        % remove the final part of empirical covariance
        if eCovF(2) < 0
            tauGrid = tauGrid_poly(1:round(size(tauGrid,1)*3/4));
            eCovF   = eCovF_poly(1:round(size(eCovF,1)*3/4));
            Cecf    = Cecf_poly(1:round(size(Cecf,1)*3/4));
        else
            tauGrid = tauGrid(1:round(size(tauGrid,1)*3/4));
            eCovF   = eCovF(1:round(size(eCovF,1)*3/4));
            Cecf    = Cecf(1:round(size(Cecf,1)*3/4));
        end
        
        tau = tauGrid(idx_in:idx_end);
        Q   = Cecf(idx_in:idx_end);
        Yo  = eCovF(idx_in:idx_end);
        
        % initialize variables for polyfit
        Yo_poly = Yo;
        tau_poly = tau;
        Q_poly = Q;
        sum_out = 1;
        max_it_out = min(4, max(2, round(size(Yo_poly,1) * 0.1)));
        it_out = 1;
        max_order = 8;
        
        % fit a polynomial to the empirical covariance with outlier removal
        l_out = min(20, idx_end - sum_out); 
        while it_out < max_it_out && sum_out~=0
            order = min(max_order, round(size(Yo_poly,1) / 3));
            p_eCov = polyfit(tau_poly, Yo_poly, order);
            eCovF_smooth_tmp = polyval(p_eCov, tau_poly);
        
            % compute outliers
            if eCovF(3) > 0
                residuals_eCovF = Yo_poly(1:l_out) - eCovF_smooth_tmp(1:l_out);
            else
                residuals_eCovF = Yo_poly(2:l_out) - eCovF_smooth_tmp(2:l_out);
            end
            std_res_eCovF = std(residuals_eCovF);
            outliers_eCovF = abs(residuals_eCovF) > 1.5 * std_res_eCovF;
            sum_out = sum(outliers_eCovF);
        
            % remove outliers
            if eCovF(3) > 0
                idx_out = find(outliers_eCovF, 1);
            else
                idx_out = find(outliers_eCovF, 1);
                idx_out = idx_out + 1;
            end
            tau_poly(idx_out) = [];
            Yo_poly(idx_out) = [];
            Q_poly(idx_out) = [];
            l_out = l_out - 1;
        
            % increment iteration
            it_out = it_out + 1;
        end
        
        % remove the initial point added to help the fit
        if eCovF(3) < 0
            tauGrid(2) = [];
            eCovF(2)   = [];
            Cecf(2)    = [];
        end
        
        % evaluate the smoothed covariance in all sampling distances
        eCovF_smooth = polyval(p_eCov, tauGrid(1:idx_end));
        
        % drop last point if it deviates > 2 std from mean of previous 3
        if length(eCovF_smooth) >= 4
            prev_three = eCovF_smooth(end-3:end-1);
            mean_prev = mean(prev_three);
            std_prev = std(prev_three);
            last_point = eCovF_smooth(end);
            if abs(last_point - mean_prev) > 2 * std_prev
                eCovF_smooth = eCovF_smooth(1:end-1);
                idx_end = idx_end - 1;
            end
        end
        
        % fix to 1e-3 the first smoothed value if negative
        if eCovF_smooth(1) < 0
            eCovF_smooth(1) = 1e-3;
        end
        
        % store the non-resampled smoothed empirical covariance
        tauGrid_orig = tauGrid;
        
        % resample for parametric fits if number of observations is low
        min_obs = idx_end_min * 2;   % threshold for resampling (30 points)
        if length(Yo) < min_obs
            fine_step = tauGrid(idx_end) / (min_obs - 1);   % finer grid step
            tau_fine = 0:fine_step:tauGrid(idx_end);
            eCovF_smooth_fine = polyval(p_eCov, tau_fine)';
            Q_fine = interp1(tauGrid(1:idx_end), Cecf(1:idx_end), tau_fine, 'linear', 'extrap');
            tau = tau_fine;
            Yo = eCovF_smooth_fine;
            Q = Q_fine;
            eCovF_smooth = eCovF_smooth_fine;
            tauGrid = tauGrid(1):fine_step:tauGrid(end);
        end
        
        % identify zero-crossing point in smoothed curve
        idxZero_smt = find(eCovF_smooth < 0, 1, 'first');
        if idxZero_smt == 1
            tauZero_smt = tauGrid(2)/2;   % fallback
        else
            tauZero_smt = mean(tauGrid(idxZero_smt-1 : idxZero_smt));
        end
        
        % identify zero-crossing point in outlier-free observations
        idx_pos_poly = find(Yo_poly < 0, 2, 'first');
        if length(idx_pos_poly) < 2 || idx_pos_poly(2) - idx_pos_poly(1) <= 1
            tau_fit_b = tau(1:min(3, length(tau)));   % fallback
            cov_fit_b = Yo(1:min(3, length(tau)));
        else
            tau_fit_b = tau_poly(1:idx_pos_poly(2)-1); tau_fit_b(idx_pos_poly(1)) = [];
            cov_fit_b = Yo_poly(1:idx_pos_poly(2)-1); cov_fit_b(idx_pos_poly(1)) = [];
        end
        
        % objective functions for Matlab built-in least squares function
        fun1 = @(c) (1./sqrt(Q)) .* (Yo - mCovF1(c, tau));
        fun2 = @(c) (1./sqrt(Q)) .* (Yo - mCovF2(c, tau));
        fun3 = @(c) (1./sqrt(Q)) .* (Yo - mCovF3(c, tau));
        
        % max deviation from a-priori values
        thrs_bounds = 0.05;
        
        % standardized Cecf
        Cecf_std = 1 - (Cecf(idx_in:idx_end) - min(Cecf(idx_in:idx_end))) / (max(Cecf(idx_in:idx_end)) - min(Cecf(idx_in:idx_end)));
        
        % -- Gaussian --
        % a-priori values for non-linear LS
        % c1
        if eCovF(2) < 0 && eCovF_smooth(2) < 0
            % fallback
            if eCovF_smooth(1) > 0 && eCovF_smooth(1) < eCovF(1)
                c1_app(1,1) = eCovF_smooth(1);
            elseif eCovF_smooth(1) > 0 && eCovF_smooth(1) > eCovF(1)
                c1_app(1,1) = min(min(eCovF_smooth(1)/2, eCovF(1)), eCovF(1)/2);
            else
                c1_app(1,1) = 1e-3;
            end
        else
            % ideal case
            c1_app(1,1) = max(eCovF(2), eCovF_smooth(2));
        end
        
        % c2
        if eCovF(2) < 0
            % fallback
            diff_eCovF = diff(eCovF_smooth);
            idx_min = find(diff_eCovF(1:end-1) < 0 & diff_eCovF(2:end) > 0, 1, 'first') + 1;
            if isempty(idx_min)
                idx_min = find(eCovF_smooth < 0, 1, 'first');
            end
            tau_min = tau(idx_min);
            f = 0.5;   % decay fraction
            c1_app(2,1) = log(1/f) / (tau_min^2);
        else
            % ideal case
            % option 1: first non-negative covariances
            tau_fit = tauGrid(idx_in:idxZero_smt-1);
            cov_fit = eCovF_smooth(idx_in:idxZero_smt-1);
            % ensure cov_fit is positive to avoid complex log
            valid_idx = cov_fit > 0;
            % need at least 2 points for polyfit
            if sum(valid_idx) >= 2   
                x_c2 = tau_fit(valid_idx).^2;
                y_c2 = log(cov_fit(valid_idx));
                p_c2 = polyfit(x_c2, y_c2, 1);
                p_c2a = -p_c2(1);
            else
                % mark as invalid
                p_c2a = NaN;   
            end
        
            % option 2: fit a straight line to get the slope
            valid_idx_b = cov_fit_b > 0;
            if sum(valid_idx_b) >= 2
                x_c2b = tau_fit_b(valid_idx_b).^2;
                y_c2b = log(cov_fit_b(valid_idx_b));
                p_corrL = polyfit(x_c2b, y_c2b, 1);
                p_c2b = -p_corrL(1);
            else
                % mark as invalid
                p_c2b = NaN;   
            end
        
            % select the minimum positive, real slope
            if isreal(p_c2a) && p_c2a > 0 && isreal(p_c2b) && p_c2b > 0
                c1_app(2,1) = min(p_c2a, p_c2b);
            elseif isreal(p_c2a) && p_c2a > 0
                c1_app(2,1) = p_c2a;
            elseif isreal(p_c2b) && p_c2b > 0
                c1_app(2,1) = p_c2b;
            else
                % fallback: use the same logic as eCovF(2) < 0
                diff_eCovF = diff(eCovF_smooth);
                idx_min = find(diff_eCovF(1:end-1) < 0 & diff_eCovF(2:end) > 0, 1, 'first') + 1;
                if isempty(idx_min)
                    idx_min = find(eCovF_smooth < 0, 1, 'first');
                end
                if ~isempty(idx_min)
                    tau_min = tau(idx_min);
                    f = 0.5;
                    c1_app(2,1) = log(1/f) / (tau_min^2);
                else
                    % ultimate fallback: use a small positive value
                    c1_app(2,1) = 1e-3;
                end
            end
        end
        
        % lower and upper bound for lsqnonlin
        lb_1 = c1_app - [1, 2] .* thrs_bounds .* c1_app;
        lb_1(lb_1 < 0) = 0;
        ub_1 = c1_app + [1, 2] .* thrs_bounds .* c1_app;
        
        % non-linear LS
        [c1, ~] = lsqnonlin(fun1, c1_app, lb_1, ub_1);
        
        % compute the score
        scoreCov_raw(1) = sum(((mCovF1(c1, tauGrid_orig(idx_in:idx_end)) - eCovF(idx_in:idx_end)) .* Cecf_std).^2);
        scoreCov_smt(1) = sum(((mCovF1(c1, tauGrid_orig(idx_in:idx_end)) - eCovF_smooth(1:idx_end-1)) .* Cecf_std).^2);
        
        % -- Gaussian with bell --
        % a-priori values for non-linear LS
        % c1
        if (eCovF(2) < 0 && eCovF_smooth(2) < 0) && eCovF_smooth(1) < eCovF(1)
            % fallback
            if eCovF_smooth(1) > 0
                c2_app(1,1) = eCovF_smooth(1);
            else
                c2_app(1,1) = 1e-3;
            end
        elseif (eCovF(2) < 0 && eCovF_smooth(2) < 0) && eCovF_smooth(1) > eCovF(1)
            c2_app(1,1) = min(min(eCovF_smooth(1)/2, eCovF(1)), eCovF(1)/2);
        else
            % ideal case
            c2_app(1,1) = max(max(eCovF(2), eCovF_smooth(2)), 1e-3);
        end
        
        % c2
        if eCovF(2) < 0
            % fallback
            diff_eCovF = diff(eCovF_smooth);
            idx_min = find(diff_eCovF(1:end-1) < 0 & diff_eCovF(2:end) > 0, 1, 'first') + 1;
            if isempty(idx_min)
                idx_min = find(eCovF_smooth < 0, 1, 'first');
            end
            tau_min = tau(idx_min);
            f = 0.5;
            c2_app(2,1) = log(1/f) / (tau_min^2);
        else
            % ideal case
            idxBellStart = idxZero_smt;
            idxBellEnd = find(eCovF_smooth(idxBellStart:end) > 0, 1, 'first') + idxBellStart - 1;
            if isempty(idxBellStart), idxBellStart = idxZero; end
            if isempty(idxBellEnd), idxBellEnd = length(tauGrid); end
            if idxBellEnd > length(tauGrid), idxBellEnd = length(tau); end
            tau_fit_bell_2 = [tau(1), tau(idxBellEnd)];
            cov_fit_bell_2 = [eCovF_smooth(1), eCovF_smooth(idxBellEnd)];
            if all(cov_fit_bell_2 > 0) && length(cov_fit_bell_2) >= 2
                x_bell_c2 = tau_fit_bell_2.^2;
                y_bell_c2 = log(cov_fit_bell_2);
                p_bell_c2 = polyfit(x_bell_c2, y_bell_c2, 1);
                c2_app(2,1) = -p_bell_c2(1);
                if ~isreal(c2_app(2,1)) || c2_app(2,1) <= 0
                    % fallback
                    c2_app(2,1) = log(1/0.5) / (tauGrid(idxBellEnd)^2);   
                end
            else
                % fallback
                c2_app(2,1) = log(1/0.5) / (tauGrid(idxBellEnd)^2);   
            end
        end
        
        % c3
        if eCovF(2) < 0
            % fallback
            gaussian_term = c2_app(1) * exp(-c2_app(2) * tau_min^2);
            if gaussian_term > 0
                c2_app(3) = (1/tau_min^2) * (1 - eCovF_smooth(idx_min)/gaussian_term);
                % ensure non-negative for bell shape
                c2_app(3) = max(c2_app(3), 0);   
            else
                % default to no bell
                c2_app(3) = 0;   
            end
        else
            % ideal case
            tau_fit_bell = tauGrid(idxBellStart:idxBellEnd);
            cov_fit_bell = eCovF_smooth(idxBellStart:idxBellEnd);
            % need at least 3 points for quadratic fit
            if sum(cov_fit_bell > 0) >= 3   
                p_bell = polyfit(tau_fit_bell.^2, cov_fit_bell, 2);
                % avoid division by zero
                if abs(p_bell(2)) > 1e-6   
                    c2_app(3,1) = -2 * p_bell(1) / p_bell(2);
                    % ensure non-negative
                    c2_app(3,1) = max(c2_app(3,1), 0);   
                else
                    % default to no bell
                    c2_app(3,1) = 0;   
                end
            else
                c2_app(3,1) = 0;   % default to no bell
            end
        end
        
        % lower and upper bound for lsqnonlin
        if eCovF(2) < 0
            % fallback
            lb_2 = [c2_app(1) - thrs_bounds .* c2_app(1), c2_app(2) - thrs_bounds .* c2_app(2), c2_app(3) - thrs_bounds .* c2_app(3)];
            ub_2 = [c2_app(1) + thrs_bounds .* c2_app(1), c2_app(2) + thrs_bounds .* c2_app(2), c2_app(3) + thrs_bounds .* c2_app(3)];
        else
            % ideal case
            lb_2 = [c2_app(1) - thrs_bounds .* c2_app(1), 0, 0];
            ub_2 = [c2_app(1) + thrs_bounds .* c2_app(1), Inf, Inf];
        end
        
        % non-linear LS
        [c2, ~] = lsqnonlin(fun2, c2_app, lb_2, ub_2);
        
        % compute the score
        scoreCov_raw(2) = sum(((mCovF2(c2, tauGrid_orig(idx_in:idx_end)) - eCovF(idx_in:idx_end)) .* Cecf_std).^2);
        scoreCov_smt(2) = sum(((mCovF2(c2, tauGrid_orig(idx_in:idx_end)) - eCovF_smooth(1:idx_end-1)) .* Cecf_std).^2);
        
        % -- Gaussian with cosine --
        % a-priori values for non-linear LS
        c3_app(1,1) = c2_app(1,1);
        c3_app(2,1) = c2_app(2,1);
        
        % c3
        if eCovF(2) < 0
            % fallback 
            diff_eCovF = diff(eCovF_smooth);
            idx_peaks = find(diff_eCovF(1:end-1) > 0 & diff_eCovF(2:end) < 0) + 1;
            if length(idx_peaks) < 2
                idx_zeros = find(eCovF_smooth(1:end-1) .* eCovF_smooth(2:end) < 0);
                if length(idx_zeros) >= 2
                    half_periods = diff(tau(idx_zeros));
                    T = 2 * mean(half_periods);
                else
                    T = tau(end) / 2;
                end
            else
                periods = diff(tau(idx_peaks));
                T = mean(periods);
            end
            if T > 0
                c3_app(3) = 4 / T;
            else
                % default period
                c3_app(3) = 1 / (dtau * 2); 
            end
        else
            % ideal case
            c3_app(3,1) = 1 / tauZero_smt;
        end
        
        % lower and upper bound for lsqnonlin
        if eCovF(2) < 0
            % fallback
            lb_3 = [c3_app(1) - thrs_bounds .* c3_app(1), c3_app(2) - thrs_bounds .* c3_app(2), c3_app(3) - thrs_bounds .* c3_app(3)];
            ub_3 = [c3_app(1) + thrs_bounds .* c3_app(1), c3_app(2) + thrs_bounds .* c3_app(2), c3_app(3) + thrs_bounds .* c3_app(3)];
        else
            % ideal case
            lb_3 = [c3_app(1) - thrs_bounds .* c3_app(1), 0, 0]; 
            ub_3 = [c3_app(1) + thrs_bounds .* c3_app(1), Inf, Inf];
        end
        
        % non-linear LS
        [c3, ~] = lsqnonlin(fun3, c3_app, lb_3, ub_3);
        
        % compute the score
        scoreCov_raw(3) = sum(((mCovF3(c3, tauGrid_orig(idx_in:idx_end)) - eCovF(idx_in:idx_end)) .* Cecf_std).^2);
        scoreCov_smt(3) = sum(((mCovF3(c3, tauGrid_orig(idx_in:idx_end)) - eCovF_smooth(1:idx_end-1)) .* Cecf_std).^2);
        
        % get total score
        scoreCov_tot = scoreCov_raw + scoreCov_smt;
        
        % get min values and indices
        [min_mCovF_raw, idx_mCovF_raw] = min(scoreCov_raw);
        [~, idx_mCovF_smt] = min(scoreCov_smt);
        [~, idx_mCovF_tot] = min(scoreCov_tot);
        
        % choice of the empirical model
        possib = [1 2 3];
        possib(idx_mCovF_raw) = [];
        if 2 * min_mCovF_raw < min(scoreCov_raw(possib))
            idx_mCovF = idx_mCovF_raw;
        elseif (idx_mCovF_raw ~=1 && idx_mCovF_smt ~=1) && idx_mCovF_tot ~=1
            idx_mCovF = idx_mCovF_tot;
        else
            % max correlation length between gau and gau + bell
            zeroF1 = find(mCovF1(c1, tauGrid) < 1e-6, 1, 'first');
            zeroF2a = find(mCovF2(c2, tauGrid) < 1e-6, 1, 'first');
            zeroF2b = find(mCovF2(c2, tauGrid(zeroF2a:end)) > -1e-6, 1, 'first') + zeroF2a - 1;
            if zeroF1 > zeroF2b
                idx_mCovF = 1;
            else
                idx_mCovF = 2;
            end
        end
        
        % get corresponding covariance model
        if idx_mCovF == 1
            mCovFin = mCovF1;
            cFin = c1;
        elseif idx_mCovF == 2
            mCovFin = mCovF2;
            cFin = c2;
        elseif idx_mCovF == 3
            mCovFin = mCovF3;
            cFin = c3;
        end
        
        % get the variance of the noise
        varNoise_cov = eCovF(1) - mCovFin(cFin, 0);
        
        % update varNoise for next iteration
        if varNoise_cov > 0 && ~isnan(varNoise_cov)
            varNoise = varNoise_cov;
        else
            warning('Invalid varNoise_cov (%.2f) for PS %d, iteration %d. Retaining previous noise variance: %.2f mm^2.', ...
                varNoise_cov, PSidIN_AOI(i), ii, varNoise);
            varNoise = varNoise_initial;
        end
        
        % do not apply collocation if the amplitude of the signal is way
        % smaller than the amplitude of the noise and correlation is short
        % normalize covariance to get autocorrelation
        rho = mCovFin(cFin, tauGrid) / mCovFin(cFin, 0);
        
        % find decorrelation time (first lag where rho drops below 1/e)
        tau_c_idx = find(rho < exp(-1), 1, 'first');
        if isempty(tau_c_idx)
            tau_c = length(tauGrid);
        else
            tau_c = tauGrid(tau_c_idx);
        end
        
        % signal-to-noise ratio
        NigNoiR = mCovFin(cFin, 0) / varNoise_cov;
        
        % threshold for tau_c
        tau_threshold = 2 * (1 / sqrt(NigNoiR)) * dtau;
        
        % condition
        if NigNoiR < 0.10 && tau_c < tau_threshold
            idx_mCovF = 4;
            mCovFin = @(c, tau) zeros(size(tau));
            cFin = 0;
        end

        % populate output matrices
        obs_p4{i, 1} = mCovFin;
        obs_p4{i, 2} = cFin;
        obs_p4{i, 3} = varNoise_cov;
    
    
        % LOOP CONDITION
        if ii > 1
    
            if varNoise_cov < varNoise_cov_prec + 0.05 * varNoise_cov_prec && ...
                varNoise_cov > varNoise_cov_prec - 0.05 * varNoise_cov_prec
                break
            else
                varNoise_cov_prec = varNoise_cov;
            end
        
        else
            varNoise_cov_prec = varNoise_cov;
        end
    
        ii = ii + 1;
        sprintf('-- Iteration %i', ii)

    end

    % Warning for max_iterations reached
    if ii > max_iterations
        warning('Max iterations reached for PS %d without full convergence. Final noise variance: %.2f mm^2.', ...
        PSidIN_AOI(i), varNoise);
    end


    % -- 5) Collocation
    if idx_mCovF == 4
        % collocation is NOT applied
        switch coll_proc
            case 'filtering'
            % populate output matrix
            t_obs = obs_p2{i, 4}';
            final_t_rel = t_obs';
            final_t_date = t_obs + t0_date;
            coll_signal = obs_p3{i, 4};   % only splines
            coll_std = zeros(size(coll_signal));

            case 'prediction'
            % resample the splines solution
            t_est = coll_t_est';
            num_spl = obs_p3{i, 1};
            lambda = obs_p3{i, 2};
            n_params = num_spl + 2;
            jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, num_spl, ...
                        t_in, t_fin, lambda, num_sig, gS_input_path, gS_output_path, gS_job_path);
        
            % run geoSplinter_analysis with the job file
            job_execution = sprintf('%s < %s', fullfile(gS_dir, 'geoSplinter_analysis'), fullfile('.', gS_job_path, [file_out, '.job']));
            status = system(job_execution);
            if status ~= 0
                error('Error executing geoSplinter_analysis for file: %s', file_out);
            end

            % import normal matrix and normal vector
            mat_file = fullfile('.', gS_output_path, strcat(file_out, '.mat.txt'));
            fid = fopen(mat_file, 'r');
            if fid == -1
                error('Cannot open file: %s', mat_file);
            end
            try
                normal_matrix = zeros(n_params, n_params);
                normal_vector = zeros(n_params, 1);
                while ~feof(fid)
                    line = fgetl(fid);
                    if ~ischar(line) || isempty(strtrim(line))
                        continue;
                    end
                    if contains(line, 'Normal Matrix') && ~contains(line, 'as stored') && ~contains(line, 'Diagonal')
                        r = 1;
                        while r <= n_params
                            line = fgetl(fid);
                            if ~ischar(line) || isempty(strtrim(line))
                                continue;
                            end
                            values = sscanf(line, repmat('%f ', 1, n_params));
                            if length(values) == n_params
                                normal_matrix(r, :) = values';
                                r = r + 1;
                            end
                        end
                    end
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
                end
            catch err
                fclose(fid);
                error('Error reading file %s: %s', mat_file, err.message);
            end
            fclose(fid);
            
            % import knots
            hdr_file = fullfile('.', gS_output_path, strcat(file_out, '.hdr.txt'));
            fid = fopen(hdr_file, 'r');
            if fid == -1
                error('Cannot open file: %s', hdr_file);
            end
            try
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
            catch err
                fclose(fid);
                error('Error reading file %s: %s', hdr_file, err.message);
            end
            fclose(fid);
            knots = xMin:deltaX:xMax;
            if isempty(knots) || length(knots) < 2
                error('Invalid knot vector for PS %d: knots must have at least 2 elements.', PSidIN_AOI(i));
            end
            
            % apply regularization for master epoch constraint
            if ~isempty(master_idx)
                B = compute_B_spline_row(master_t_rel, knots, deltaX, n_params);
                if sum(B) == 0
                    warning('B-spline row is all zeros for master epoch for PS %d.', PSidIN_AOI(i));
                end
                lambdaR = 1e1;
                R = B' * B;
                N_modified = normal_matrix + lambdaR * R;
                x_constrained = N_modified \ normal_vector;
            else
                % fallback to unconstrained
                warning('Master epoch not found in time vector for PS %d.', PSidIN_AOI(i));
                x_constrained = normal_matrix \ normal_vector; 
            end

            % replicate geoSplinter_analysis output files for synthesis
            % write .par.txt with constrained coefficients
            constrained_par_file = fullfile(gS_output_path, strcat(file_out, '.par.txt'));
            fid = fopen(constrained_par_file, 'w');
            if fid == -1
                error('Cannot create file: %s', constrained_par_file);
            end
            try
                for r = 1:n_params
                    fprintf(fid, '%+e\n', x_constrained(r));
                end
            catch err
                fclose(fid);
                error('Error writing file %s: %s', constrained_par_file, err.message);
            end
            fclose(fid);
    
            % synthetize splines solution on t_est
            gS_est = 't_est';
            writematrix(t_est, strcat(gS_input_path, filesep, gS_est, '.txt'), 'Delimiter', 'space')
            
            % job file creation
            file_spl = file_out;                                     % output file from geoSplinter analysis
            file_est = strcat(gS_est, '.txt');                       % query times filename (with extension)
            file_syn = strcat(gS_filename, '_', ext, '_', 'est');    % output filename (without extension)
            
            jobFile_synthesis(data_dim, type_spl, file_spl, file_est, file_syn, gS_input_path, gS_output_path, gS_job_path, gS_synth_path);
            
            % run geoSplinter_synthesis with the job file
            job_execution_syn = sprintf('%s < %s', fullfile(gS_dir, 'geoSplinter_synthesis'), fullfile('.', gS_job_path, [file_syn, '.job']));
            status = system(job_execution_syn);
            if status ~= 0
                error('Error executing geoSplinter_synthesis for file: %s', file_syn);
            end
                
            % results import
            spl_est = readmatrix(fullfile(gS_synth_path, file_syn));
            spl_est(:,1) = [];

            % populate output matrix
            final_t_rel = t_est';
            final_t_date = (t_est + t0_date)';
            coll_signal = spl_est';
            coll_std = zeros(size(coll_signal));
        end

    else
        % collocation is applied
        switch coll_proc
            case 'filtering'
            % FILTERING
            % observations
            t_obs = obs_p2{i, 4}';
            n_obs = length(t_obs);
            
            % estimations
            t_est = t_obs;
            n_est = n_obs;
            
            % verify spline residual at master epoch
            if ~isempty(master_idx)
                master_residual = obs_p3{i, 5}(master_idx);
                if abs(master_residual) > 1e-1
                    warning('Spline residual is not zero at master epoch for PS %d: %.10f', PSidIN_AOI(i), master_residual);
                end
            end
            
            % covariance matrix est-obs
            T_est = repmat(t_est, 1, n_obs);
            T_obs = repmat(t_obs', n_est, 1);
            C_est_obs = obs_p4{i, 1}(obs_p4{i, 2}, abs(T_est - T_obs));
            
            % covariance matrix obs
            T_obs = repmat(t_obs', n_obs, 1);
            C_obs_obs = obs_p4{i, 1}(obs_p4{i, 2}, abs(T_obs' - T_obs));
    
            % check for positive definite covariance matrix (for inversion)
            min_eig = min(eig(C_obs_obs));
            if min_eig <= 0
                if ~exist(fullfile('MatlabFunctions', 'nearestSPD'), 'file')
                    warning('nearestSPD not found. Using diagonal adjustment for positive definiteness.');
                    C_obs_obs = C_obs_obs + 1e-6 * eye(size(C_obs_obs)); % small diagonal adjustment
                else
                    C_obs_obs = nearestSPD(C_obs_obs);
                end
            end

            % modify noise covariance matrix (Cvv)
            Cvv = obs_p4{i, 3} * eye(n_obs);
            if ~isempty(master_idx)
                master_variance = 1e-1; 
                Cvv(master_idx, master_idx) = master_variance;
            end
            
            % compute weighted mean
            weights = 1 ./ diag(Cvv);
            weighted_mean = (weights' * obs_p3{i, 5}') / sum(weights);

            % collocation solution
            coll_est = C_est_obs * (C_obs_obs + Cvv)^(-1) * (obs_p3{i, 5}' - weighted_mean);
            coll_est = coll_est + weighted_mean;
    
            case 'prediction'
            % PREDICTION
            % observations
            t_obs = obs_p2{i, 4}';
            n_obs = length(t_obs);
            
            % estimations
            t_est = coll_t_est';
            n_est = length(t_est);
            
            % covariance matrix est-obs
            T_est = repmat(t_est, 1, n_obs);
            T_obs = repmat(t_obs', n_est, 1);
            C_est_obs = obs_p4{i, 1}(obs_p4{i, 2}, abs(T_est - T_obs));
            
            % covariance matrix obs
            T_obs = repmat(t_obs', n_obs, 1);
            C_obs_obs = obs_p4{i, 1}(obs_p4{i, 2}, abs(T_obs' - T_obs));
    
            % check for positive definite covariance matrix (for inversion)
            min_eig = min(eig(C_obs_obs));
            if min_eig <= 0
                if ~exist(fullfile('MatlabFunctions', 'nearestSPD'), 'file')
                    warning('nearestSPD not found. Using diagonal adjustment for positive definiteness.');
                    C_obs_obs = C_obs_obs + 1e-6 * eye(size(C_obs_obs)); % small diagonal adjustment
                else
                    C_obs_obs = nearestSPD(C_obs_obs);
                end
            end

            % modify noise covariance matrix (Cvv)
            Cvv = obs_p4{i, 3} * eye(n_obs);
            if ~isempty(master_idx)
                master_variance = 1e-1; 
                Cvv(master_idx, master_idx) = master_variance;
            end

            % compute weighted mean
            weights = 1 ./ diag(Cvv);
            weighted_mean = (weights' * obs_p3{i, 5}') / sum(weights);
            
            % collocation solution
            coll_est = C_est_obs * (C_obs_obs + Cvv)^(-1) * (obs_p3{i, 5}' - weighted_mean);
            coll_est = coll_est + weighted_mean;
    
        end
    
    
        % -- 6) Prediction error
        coll_var = max(obs_p4{i, 1}(obs_p4{i, 2}, 0) - diag(C_est_obs * (C_obs_obs + Cvv)^(-1) * C_est_obs'), 10^(-6));
        coll_std = sqrt(coll_var');
    
        % final solution
        switch coll_proc
            case 'filtering'
            % FILTERING
            final_t_rel = obs_p2{i, 4};
            final_t_date = obs_p3{i, 3};
            coll_signal = obs_p3{i, 4} + coll_est';
    
            case 'prediction'
            % PREDICTION
            % re-compute the splines solution with the chosen number of splines
            num_spl = obs_p3{i, 1};
            lambda = obs_p3{i, 2};
            n_params = num_spl + 2;
            jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, num_spl, ...
                        t_in, t_fin, lambda, num_sig, gS_input_path, gS_output_path, gS_job_path);
        
            % run geoSplinter_analysis with the job file
            job_execution = sprintf('%s < %s', fullfile(gS_dir, 'geoSplinter_analysis'), fullfile('.', gS_job_path, [file_out, '.job']));
            status = system(job_execution);
            if status ~= 0
                error('Error executing geoSplinter_analysis for file: %s', file_out);
            end

            % import normal matrix and normal vector
            mat_file = fullfile('.', gS_output_path, strcat(file_out, '.mat.txt'));
            fid = fopen(mat_file, 'r');
            if fid == -1
                error('Cannot open file: %s', mat_file);
            end
            try
                normal_matrix = zeros(n_params, n_params);
                normal_vector = zeros(n_params, 1);
                while ~feof(fid)
                    line = fgetl(fid);
                    if ~ischar(line) || isempty(strtrim(line))
                        continue;
                    end
                    if contains(line, 'Normal Matrix') && ~contains(line, 'as stored') && ~contains(line, 'Diagonal')
                        r = 1;
                        while r <= n_params
                            line = fgetl(fid);
                            if ~ischar(line) || isempty(strtrim(line))
                                continue;
                            end
                            values = sscanf(line, repmat('%f ', 1, n_params));
                            if length(values) == n_params
                                normal_matrix(r, :) = values';
                                r = r + 1;
                            end
                        end
                    end
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
                end
            catch err
                fclose(fid);
                error('Error reading file %s: %s', mat_file, err.message);
            end
            fclose(fid);
            
            % import knots
            hdr_file = fullfile('.', gS_output_path, strcat(file_out, '.hdr.txt'));
            fid = fopen(hdr_file, 'r');
            if fid == -1
                error('Cannot open file: %s', hdr_file);
            end
            try
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
            catch err
                fclose(fid);
                error('Error reading file %s: %s', hdr_file, err.message);
            end
            fclose(fid);
            knots = xMin:deltaX:xMax;
            if isempty(knots) || length(knots) < 2
                error('Invalid knot vector for PS %d: knots must have at least 2 elements.', PSidIN_AOI(i));
            end
            
            % apply regularization for master epoch constraint
            if ~isempty(master_idx)
                B = compute_B_spline_row(master_t_rel, knots, deltaX, n_params);
                if sum(B) == 0
                    warning('B-spline row is all zeros for master epoch for PS %d.', PSidIN_AOI(i));
                end
                lambdaR = 1e1;
                R = B' * B;
                N_modified = normal_matrix + lambdaR * R;
                x_constrained = N_modified \ normal_vector;
            else
                % fallback to unconstrained
                warning('Master epoch not found in time vector for PS %d.', PSidIN_AOI(i));
                x_constrained = normal_matrix \ normal_vector; 
            end

            % replicate geoSplinter_analysis output files for synthesis
            % write .par.txt with constrained coefficients
            constrained_par_file = fullfile(gS_output_path, strcat(file_out, '.par.txt'));
            fid = fopen(constrained_par_file, 'w');
            if fid == -1
                error('Cannot create file: %s', constrained_par_file);
            end
            try
                for r = 1:n_params
                    fprintf(fid, '%+e\n', x_constrained(r));
                end
            catch err
                fclose(fid);
                error('Error writing file %s: %s', constrained_par_file, err.message);
            end
            fclose(fid);
    
            % synthetize splines solution on t_est
            gS_est = 't_est';
            writematrix(t_est, strcat(gS_input_path, filesep, gS_est, '.txt'), 'Delimiter', 'space')
            
            % job file creation
            file_spl = file_out;                                     % output file from geoSplinter analysis
            file_est = strcat(gS_est, '.txt');                       % query times filename (with extension)
            file_syn = strcat(gS_filename, '_', ext, '_', 'est');    % output filename (without extension)
            
            jobFile_synthesis(data_dim, type_spl, file_spl, file_est, file_syn, gS_input_path, gS_output_path, gS_job_path, gS_synth_path);
            
            % run geoSplinter_synthesis with the job file
            job_execution_syn = sprintf('%s < %s', fullfile(gS_dir, 'geoSplinter_synthesis'), fullfile('.', gS_job_path, [file_syn, '.job']));
            status = system(job_execution_syn);
            if status ~= 0
                error('Error executing geoSplinter_synthesis for file: %s', file_syn);
            end
                
            % results import
            spl_est = readmatrix(fullfile(gS_synth_path, file_syn));
            spl_est(:,1) = [];
    
            % populate output matrix
            final_t_rel = t_est';
            final_t_date = (t_est + t0_date)';
            coll_signal = spl_est' + coll_est';
    
        end

    end

    % RECONTRUCT THE SIGNAL
    % trend + periodicity + collocation
    switch coll_proc
        case 'filtering'
            trend_signal = obs_p1{i, 4};
            periodic_signal = obs_p2{i, 7};
        final_signal = trend_signal + periodic_signal + coll_signal;
        
        case 'prediction'
            % trend signal
            coeffs = obs_p1{i, 6};
            model_type = obs_p1{i, 7};
            switch model_type
                case 'no_model'
                    trend_signal = zeros(size(final_t_rel));
                case 'constant'
                    trend_signal = coeffs(1) * ones(size(final_t_rel));
                case 'linear'
                    trend_signal = coeffs(1) + coeffs(2) * final_t_rel;
                case 'quadratic'
                    trend_signal = coeffs(1) + coeffs(2) * final_t_rel + coeffs(3) * final_t_rel.^2;
            end
            % periodic signal
            fundamental_freq = obs_p2{i, 2};
            coeffs = obs_p2{i, 6};
            if ~isempty(fundamental_freq) && ~isempty(coeffs)
                n_freq = length(fundamental_freq);
                periodic_signal = zeros(size(final_t_rel));
                % constant term
                periodic_signal = periodic_signal + coeffs(1);
                % cosine and sine terms
                for k = 1:n_freq
                    periodic_signal = periodic_signal + ...
                        coeffs(2*k) * cos(2 * pi * fundamental_freq(k) * final_t_rel) + ...
                        coeffs(2*k+1) * sin(2 * pi * fundamental_freq(k) * final_t_rel);
                end
            else
                % no periodicity detected, set periodic signal to zero
                periodic_signal = zeros(size(final_t_rel));
            end
        % final signal
        final_signal = trend_signal + periodic_signal + coll_signal;
    end

    % DEFINE THE UNCERTAINTY OF THE ESTIMATE
    % - trend signal
    % retrieve polynomial coefficients and their standard deviations
    std_poly_final = obs_p1{i, 5};
    model_type = obs_p1{i, 7};
    t_trend = final_t_rel;
    trend_std = zeros(1, length(t_trend));
    
    switch model_type
        case 'no_model'
            % no model, zero uncertainty
            trend_std = zeros(1, length(t_trend)); 
        case 'constant'
            % only a0, constant uncertainty
            trend_std = std_poly_final(1) * ones(1, length(t_trend)); 
        case {'linear', 'quadratic'}
            % design matrix for polynomial
            if strcmp(model_type, 'linear')
                A_trend = [ones(length(t_trend), 1), t_trend'];   % [1, t]
                Cxx_trend = diag(std_poly_final(1:2).^2);   % covariance for [a0, a1]
            else % quadratic
                A_trend = [ones(length(t_trend), 1), t_trend', t_trend'.^2];   % [1, t, t^2]
                Cxx_trend = diag(std_poly_final.^2);   % covariance for [a0, a1, a2]
            end
            % propagate variance: var(y) = A * Cxx * A'
            trend_var = zeros(1, length(t_trend));
            for j = 1:length(t_trend)
                trend_var(1,j) = A_trend(j,:) * Cxx_trend * A_trend(j,:)';
            end
            trend_std = sqrt(trend_var);
    end

    % - periodic signal
    fundamental_freq = obs_p2{i, 2};
    coeffs_periodic = obs_p2{i, 6};
    t_periodic = final_t_rel;
    
    if ~isempty(fundamental_freq) && ~isempty(coeffs_periodic)
        n_freq = length(fundamental_freq);
        % design matrix for periodic signal
        A_periodic = zeros(length(t_periodic), 2 * n_freq + 1);
        A_periodic(:, 1) = 1; % constant term
        for k = 1:n_freq
            A_periodic(:, 2*k) = cos(2 * pi * fundamental_freq(k) * t_periodic);
            A_periodic(:, 2*k+1) = sin(2 * pi * fundamental_freq(k) * t_periodic);
        end
        
        % compute covariance of coefficients
        % residuals from periodic fit
        residuals_periodic = obs_p1{i, 3} - obs_p2{i, 7};   % residuals at original times
        s02_periodic = var(residuals_periodic);   % variance of residuals
        W_periodic = eye(length(obs_p2{i, 4}));
        master_idx_reg = find(abs(obs_p2{i, 4} - master_t_rel) < eps, 1);
        if ~isempty(master_idx_reg)
            master_variance = 1e-1;
            W_periodic(master_idx_reg, master_idx_reg) = 1 / master_variance;
        end
        % original design matrix at observation times
        A_periodic_orig = zeros(length(obs_p2{i, 4}), 2 * n_freq + 1);
        A_periodic_orig(:, 1) = 1;
        for k = 1:n_freq
            A_periodic_orig(:, 2*k) = cos(2 * pi * fundamental_freq(k) * obs_p2{i, 4}');
            A_periodic_orig(:, 2*k+1) = sin(2 * pi * fundamental_freq(k) * obs_p2{i, 4}');
        end
        Cxx_periodic = s02_periodic * ((A_periodic_orig' * W_periodic * A_periodic_orig) \ eye(size(A_periodic_orig, 2)));
        
        % propagate variance
        periodic_var = zeros(size(t_periodic));
        for j = 1:length(t_periodic)
            periodic_var(1,j) = A_periodic(j,:) * Cxx_periodic * A_periodic(j,:)';
        end
        periodic_std = sqrt(periodic_var);
    else
        % no periodic signal, zero uncertainty
        periodic_std = zeros(1, length(t_periodic)); 
    end

    % - spline + collocation signal
    t_spline = final_t_rel;
    n_params = num_spl + 2;
    knots = knots_csn;
    deltaX = abs(knots(3) - knots(2));

    % compute spline design matrix
    A_spline = zeros(length(t_spline), n_params);
    for n = 1:length(t_spline)
        t = t_spline(n);
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
        alpha(1) = phi_4(csi_i + 1);
        alpha(2) = phi_3(csi_i);
        alpha(3) = phi_3(1 - csi_i);
        alpha(4) = phi_4(2 - csi_i);
        for s = -1:2
            u = i_i + s + 1;
            if u >= 1 && u <= n_params
                A_spline(n, u) = alpha(s + 2);
            end
        end
    end

    % covariance matrix of spline coefficients assuming uncorrelated coeffs
    Cxx_spline = diag(std_splParms_full.^2);
    spline_var = zeros(1, length(t_spline));
    % propagate variance
    for j = 1:length(t_spline)
        spline_var(1,j) = A_spline(j,:) * Cxx_spline * A_spline(j,:)';
    end
    spline_std = sqrt(spline_var);

    % sum the two components
    tot_coll_std = sqrt(spline_std.^2 + coll_std.^2);

    % - final variance and std (Conditional Sequential Logic)
        % check for positiveness
        trend_std = max(trend_std, 0);
        periodic_std = max(periodic_std, 0);
        tot_coll_std = max(tot_coll_std, 0);
    if sum(coll_std) == 0 || all(isnan(coll_std))
        % collocation was skipped
        final_std = spline_std;
        out_coll_std = spline_std; 
    else
        % collocation was performed
        final_std = coll_std;
        out_coll_std = coll_std;
    end


    % populate output matrices
    obs_p5{i, 1} = final_t_rel;
    obs_p5{i, 2} = final_t_date;
    obs_p5{i, 3} = [trend_signal; trend_std];
    obs_p5{i, 4} = [periodic_signal; periodic_std];
    obs_p5{i, 5} = [coll_signal; out_coll_std];
    obs_p5{i, 6} = [final_signal; final_std];


    % clean up all temporary files for this PS point
    delete(fullfile(gS_output_path, [file_out, '*']));
    delete(fullfile(gS_job_path, [file_out, '.job']));
    delete(fullfile(gS_synth_path, [gS_filename, '_', ext, '_est*']));
    fclose('all');


    disp('Processing completed.')

end




%%           -------- END OF SCRIPT -------
disp('===========================================');
disp('--- PS time series modelling started... ---');
disp('===========================================');