function [coeffs, model_func] = fit_exp_auto(lags, cov_vals, variances)

% fit_exp_auto Fits an Exponential model using polynomial smoothing and outlier removal
%
% Inputs:
%   lags     - vector of lag distances (tauGrid)
%   cov_vals - vector of covariance values (eCovF)
%   variances- vector of estimation variances/cofactors (Cecf)

% 1. Input Management
% Ensure column vectors
tauGrid = lags(:);
eCovF   = cov_vals(:);

if nargin < 3 || isempty(variances)
    Cecf = ones(size(eCovF));
else
    Cecf = variances(:);
end

% Ensure enough data points to proceed
if length(tauGrid) < 4
    % Fallback: Simple fit if not enough points for the robust procedure
    model_func = @(c, tau) c(1) .* exp(-c(2) .* tau);
    coeffs = [max(eCovF); 1e-5]; 
    return;
end

% 2. Define Model
% Gaussian: c(1) = Amplitude, c(2) = Decay
mCovF1 = @(c, tau) c(1) .* exp(-c(2) .* tau);
model_func = mCovF1;

% 3. Pre-processing

% create temporary copies for polynomial fitting
tauGrid_poly = tauGrid;
eCovF_poly = eCovF;
Cecf_poly = Cecf;

% check if second value is negative and add constraint in temporary copies
if length(eCovF) > 2 && eCovF(2) < 0
    x = [tauGrid_poly(2), tauGrid_poly(3)];
    y = [eCovF_poly(2), eCovF_poly(3)];
    tau0_interp = interp1(x, y, 0, 'linear', 'extrap');
    tau0_interp = max(tau0_interp, 0);
    
    tauGrid_poly = [tauGrid_poly(1); 0; tauGrid_poly(2:end)];
    eCovF_poly   = [eCovF_poly(1); tau0_interp; eCovF_poly(2:end)];
    Cecf_poly    = [Cecf_poly(1); Cecf_poly(1); Cecf_poly(2:end)];
end

% find the tau coordinates for which the covariance function becomes 0
idxZero = find(eCovF < 0, 1, 'first');
if isempty(idxZero)
    % Fallback if no zero crossing found
    idxZero = length(eCovF); 
    tauZero = tauGrid(end);
else
    tauZero = mean(tauGrid(idxZero-1 : idxZero));
end

% extract covariance points to be used in interpolation
idx_in = 2; % Skip lag 0 (nugget)
idx_end_min = 15;    
idx_end = max(idx_end_min, round(size(Cecf,1)/2));
idx_end = min(idx_end, length(tauGrid)); % Safety cap

% remove the final part of empirical covariance
trim_factor = 3/4;
lim_idx = round(size(tauGrid_poly,1)*trim_factor);

if length(eCovF) > 1 && eCovF(2) < 0
    tauGrid_w = tauGrid_poly(1:min(lim_idx, end));
    eCovF_w   = eCovF_poly(1:min(lim_idx, end));
    Cecf_w    = Cecf_poly(1:min(lim_idx, end));
else
    tauGrid_w = tauGrid(1:min(lim_idx, end));
    eCovF_w   = eCovF(1:min(lim_idx, end));
    Cecf_w    = Cecf(1:min(lim_idx, end));
end

% Prepare data for polyfit
% Ensure indices are within bounds of the trimmed vectors
curr_len = length(tauGrid_w);
idx_use_end = min(idx_end, curr_len);

if idx_in > idx_use_end
     % Not enough data after trimming
     coeffs = [max(eCovF); 1e-5]; return;
end

tau = tauGrid_w(idx_in:idx_use_end);
Q   = Cecf_w(idx_in:idx_use_end);
Yo  = eCovF_w(idx_in:idx_use_end);

% initialize variables for polyfit
Yo_poly = Yo;
tau_poly = tau;
Q_poly = Q;
sum_out = 1;
max_it_out = min(4, max(2, round(size(Yo_poly,1) * 0.1)));
it_out = 1;
max_order = 8;

% 4. Polynomial Smoothing & Outlier Removal
p_eCov = [];
while it_out < max_it_out && sum_out~=0 && length(tau_poly) > max_order
    l_out = min(20, length(tau_poly)); % Safety on index
    order = min(max_order, round(size(Yo_poly,1) / 3));
    
    if order < 1, order = 1; end
    
    p_eCov = polyfit(tau_poly, Yo_poly, order);
    eCovF_smooth_tmp = polyval(p_eCov, tau_poly);

    % compute outliers
    if length(eCovF) > 2 && eCovF(3) > 0
        % Check bounds
        len_check = min(length(Yo_poly), l_out);
        residuals_eCovF = Yo_poly(1:len_check) - eCovF_smooth_tmp(1:len_check);
    else
        len_check = min(length(Yo_poly), l_out);
        if len_check > 1
            residuals_eCovF = Yo_poly(2:len_check) - eCovF_smooth_tmp(2:len_check);
        else
            residuals_eCovF = 0;
        end
    end
    
    std_res_eCovF = std(residuals_eCovF);
    if std_res_eCovF == 0, std_res_eCovF = 1; end % Avoid div by zero
    
    outliers_eCovF = abs(residuals_eCovF) > 1.5 * std_res_eCovF;
    sum_out = sum(outliers_eCovF);

    % remove outliers
    if sum_out > 0
        idx_out = find(outliers_eCovF, 1);
        tau_poly(idx_out) = [];
        Yo_poly(idx_out) = [];
        Q_poly(idx_out) = [];
    end

    it_out = it_out + 1;
end

% If polyfit failed or no data left, return default
if isempty(p_eCov)
     coeffs = [max(eCovF); 1e-5]; return;
end

% evaluate the smoothed covariance in all sampling distances
% Handle trimming of original grid
tauGrid_used = tauGrid(1:min(idx_end, length(tauGrid)));
eCovF_smooth = polyval(p_eCov, tauGrid_used);

% drop last point if it deviates > 2 std from mean of previous 3
if length(eCovF_smooth) >= 4
    prev_three = eCovF_smooth(end-3:end-1);
    mean_prev = mean(prev_three);
    std_prev = std(prev_three);
    last_point = eCovF_smooth(end);
    if abs(last_point - mean_prev) > 2 * std_prev
        eCovF_smooth = eCovF_smooth(1:end-1);
        % Update parallel arrays
        tauGrid_used = tauGrid_used(1:end-1);
    end
end

% fix to 1e-3 the first smoothed value if negative
if eCovF_smooth(1) < 0
    eCovF_smooth(1) = 1e-3;
end

% 5. Resampling (for parametric fits if obs low)
min_obs = idx_end_min * 2;   
if length(Yo) < min_obs && length(tauGrid_used) > 1
    fine_step = tauGrid_used(end) / (min_obs - 1);
    tau_fine = 0:fine_step:tauGrid_used(end);
    eCovF_smooth_fine = polyval(p_eCov, tau_fine)';
    
    % Interpolate Q on fine grid
    len_Q = min(length(tauGrid), length(Cecf));
    Q_fine = interp1(tauGrid(1:len_Q), Cecf(1:len_Q), tau_fine, 'linear', 'extrap');
    
    tau = tau_fine(:);
    Yo = eCovF_smooth_fine(:);
    Q = Q_fine(:);
    eCovF_smooth = eCovF_smooth_fine;
else
    % Use the smoothed values on the original grid
    tau = tauGrid_used;
    Yo = eCovF_smooth;
    % Q needs to match size
    Q = Cecf(1:length(Yo));
end

% Identify zero-crossing point in smoothed curve
idxZero_smt = find(eCovF_smooth < 0, 1, 'first');
if isempty(idxZero_smt)
     idxZero_smt = length(eCovF_smooth) + 1; 
end

% 6. Parameter Initialization
c1_app = zeros(2,1);
thrs_bounds = 0.05;

% --- Guess c1 (Amplitude) ---
if length(eCovF) > 1 && eCovF(2) < 0 && length(eCovF_smooth) > 1 && eCovF_smooth(2) < 0
    if eCovF_smooth(1) > 0 && eCovF_smooth(1) < eCovF(1)
        c1_app(1,1) = eCovF_smooth(1);
    elseif eCovF_smooth(1) > 0 && eCovF_smooth(1) > eCovF(1)
        c1_app(1,1) = min(min(eCovF_smooth(1)/2, eCovF(1)), eCovF(1)/2);
    else
        c1_app(1,1) = 1e-3;
    end
else
    % Ideal case
    val1 = 0; val2 = 0;
    if length(eCovF) > 1, val1 = eCovF(2); end
    if length(eCovF_smooth) > 1, val2 = eCovF_smooth(2); end
    c1_app(1,1) = max(val1, val2);
end

% --- Guess c2 (Decay) ---

% Option 1: Use log slope of smoothed data
limit_idx = min(length(eCovF_smooth), idxZero_smt-1);
if limit_idx < 2
    tau_fit = tau(1:min(end,3));
    cov_fit = Yo(1:min(end,3));
else
    tau_fit = tau(1:limit_idx);
    cov_fit = eCovF_smooth(1:limit_idx);
end

valid_idx = 1 : min([sum(cov_fit > 0), sum(cov_fit > 0.1 * eCovF_smooth(1)), 31]);

if sum(valid_idx) >= 2
    x_c2 = tau_fit(valid_idx);
    y_c2 = log(cov_fit(valid_idx));
    p_c2 = polyfit(x_c2, y_c2, 1);
    p_c2a = -p_c2(1);
else
    p_c2a = NaN;
end

% Assign c2
if isreal(p_c2a) && p_c2a > 0
    c1_app(2,1) = p_c2a;
else
    % Fallback: 1/e decay guess
    target = c1_app(1,1) * 0.37;
    [~, idx_d] = min(abs(eCovF_smooth - target));
    if idx_d == 1 && length(tau) > 1, idx_d = 2; end
    if idx_d <= length(tau)
        c1_app(2,1) = 1 / (tau(idx_d)^2);
    else
        c1_app(2,1) = 1e-5; 
    end
end

% 7. Optimization
% Objective function (Weighted Least Squares)
fun1 = @(c) (1./sqrt(Q)) .* (Yo - mCovF1(c, tau));

% Bounds (Using [1;2] to fix dimension error)
bound_factors = [1; 2]; 

lb_1 = c1_app - bound_factors .* thrs_bounds .* c1_app;
lb_1(lb_1 < 0) = 0;
ub_1 = c1_app + bound_factors .* thrs_bounds .* c1_app;

options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'trust-region-reflective');

try
    coeffs = lsqnonlin(fun1, c1_app, lb_1, ub_1, options);
catch
    coeffs = c1_app;
end

coeffs = coeffs(:);
end