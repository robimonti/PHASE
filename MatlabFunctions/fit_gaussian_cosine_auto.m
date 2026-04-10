function [coeffs, model_func] = fit_gaussian_cosine_auto(lags, cov_vals, variances)

% fit_gaussian_cosine_auto Fits a Gaussian * Cosine model following 
% polynomial smoothing and outlier removal logic.
%
% Inputs:
%   lags      - vector of lag distances (tauGrid)
%   cov_vals  - vector of covariance values (eCovF)
%   variances - vector of estimation variances/cofactors (Cecf)

% --- 1. Input Management ---
tauGrid = lags(:);
eCovF   = cov_vals(:);

if nargin < 3 || isempty(variances)
    Cecf = ones(size(eCovF));
else
    Cecf = variances(:);
end

% --- 2. Define Model ---
% Gaussian * Cosine: c(1)=Amplitude, c(2)=Decay, c(3)=Zero-crossing factor
mCovF3 = @(c, tau) c(1) .* exp(-c(2) .* tau.^2) .* cos(c(3) .* tau .* pi/2);
model_func = mCovF3;

if length(tauGrid) < 5
    coeffs = [max(eCovF); 1e-5; 1/tauGrid(end)]; 
    return;
end

% --- 3. Pre-processing & Constraints ---
tauGrid_poly = tauGrid;
eCovF_poly = eCovF;
Cecf_poly = Cecf;

% Check if second value is negative and add constraint (interpolation to zero)
if length(eCovF) > 2 && eCovF(2) < 0
    x = [tauGrid_poly(2), tauGrid_poly(3)];
    y = [eCovF_poly(2), eCovF_poly(3)];
    tau0_interp = max(interp1(x, y, 0, 'linear', 'extrap'), 0);
    
    tauGrid_poly = [tauGrid_poly(1); 0; tauGrid_poly(2:end)];
    eCovF_poly   = [eCovF_poly(1); tau0_interp; eCovF_poly(2:end)];
    Cecf_poly    = [Cecf_poly(1); Cecf_poly(1); Cecf_poly(2:end)];
end

% Trim final part (3/4 rule)
idx_end_min = 15;    
idx_end = max(idx_end_min, round(size(Cecf,1)/2));
idx_end = min(idx_end, length(tauGrid));

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

idx_in = 2; 
curr_len = length(tauGrid_w);
idx_use_end = min(idx_end, curr_len);

tau = tauGrid_w(idx_in:idx_use_end);
Q   = Cecf_w(idx_in:idx_use_end);
Yo  = eCovF_w(idx_in:idx_use_end);

% --- 4. Iterative Polynomial Smoothing & Outlier Removal ---
Yo_poly = Yo; tau_poly = tau; Q_poly = Q;
sum_out = 1; it_out = 1; max_order = 8;
max_it_out = min(4, max(2, round(size(Yo_poly,1) * 0.1)));
p_eCov = [];

while it_out < max_it_out && sum_out~=0 && length(tau_poly) > max_order
    l_out = min(20, length(tau_poly)); 
    order = min(max_order, round(size(Yo_poly,1) / 3));
    if order < 1, order = 1; end
    
    p_eCov = polyfit(tau_poly, Yo_poly, order);
    eCovF_smooth_tmp = polyval(p_eCov, tau_poly);

    % Compute outliers
    if length(eCovF) > 2 && eCovF(3) > 0
        len_check = min(length(Yo_poly), l_out);
        residuals_eCovF = Yo_poly(1:len_check) - eCovF_smooth_tmp(1:len_check);
    else
        len_check = min(length(Yo_poly), l_out);
        residuals_eCovF = (len_check > 1) * (Yo_poly(2:len_check) - eCovF_smooth_tmp(2:len_check));
    end
    
    std_res = std(residuals_eCovF);
    if std_res == 0, std_res = 1; end 
    outliers = abs(residuals_eCovF) > 1.5 * std_res;
    sum_out = sum(outliers);

    if sum_out > 0
        idx_out = find(outliers, 1);
        tau_poly(idx_out) = []; Yo_poly(idx_out) = []; Q_poly(idx_out) = [];
    end
    it_out = it_out + 1;
end

if isempty(p_eCov), coeffs = [max(eCovF); 1e-5; 1/tauGrid(end)]; return; end

% --- 5. Evaluate Smooth Curve & Resample ---
tauGrid_used = tauGrid(1:min(idx_end, length(tauGrid)));
eCovF_smooth = polyval(p_eCov, tauGrid_used);

% Drop last point if outlier relative to previous 3
if length(eCovF_smooth) >= 4
    prev_three = eCovF_smooth(end-3:end-1);
    if abs(eCovF_smooth(end) - mean(prev_three)) > 2 * std(prev_three)
        eCovF_smooth = eCovF_smooth(1:end-1);
        tauGrid_used = tauGrid_used(1:end-1);
    end
end

if eCovF_smooth(1) < 0, eCovF_smooth(1) = 1e-3; end

% Resampling if number of observations is low
min_obs = idx_end_min * 2;   
if length(Yo) < min_obs && length(tauGrid_used) > 1
    fine_step = tauGrid_used(end) / (min_obs - 1);
    tau_fine = (0:fine_step:tauGrid_used(end))';
    eCovF_smooth_fine = polyval(p_eCov, tau_fine);
    Q_fine = interp1(tauGrid(1:min(length(tauGrid), length(Cecf))), ...
                     Cecf(1:min(length(tauGrid), length(Cecf))), tau_fine, 'linear', 'extrap');
    tau = tau_fine; Yo = eCovF_smooth_fine; Q = Q_fine;
else
    tau = tauGrid_used; Yo = eCovF_smooth; Q = Cecf(1:length(Yo));
end

% --- 7. Limit Data to Second Zero-Crossing for LS ---
% Find indices where the sign changes
zero_crossings = find(diff(sign(Yo)) ~= 0);

if length(zero_crossings) >= 2
    % Cut data at the second zero-crossing to capture the full first oscillation
    idx_cut = zero_crossings(2);
    
    % Safety: Ensure we have enough points (at least 5-6 points past the min)
    [~, min_idx_local] = min(Yo(1:idx_cut));
    if idx_cut < min_idx_local + 2
         idx_cut = min(length(Yo), min_idx_local + 10);
    end

    tau_fit = tau(1:idx_cut);
    Yo_fit  = Yo(1:idx_cut);
    Q_fit   = Q(1:idx_cut);
else
    % Fallback: use all available data if oscillation is not clear
    tau_fit = tau;
    Yo_fit  = Yo;
    Q_fit   = Q;
end

% --- 8. Robust Parameter Guessing ---
c_app = zeros(3,1);
c_app(1) = max(Yo_fit(1), 1e-3); % Amplitude (c1)

% Frequency (c3): Based on first zero crossing
if ~isempty(zero_crossings)
    tauZero = tau(zero_crossings(1));
    c_app(3) = 1 / tauZero;
else
    c_app(3) = 1 / (tau_fit(end)/2); % Rough fallback
end

% Decay (c2): CRITICAL CHANGE
% Use the magnitude of the first MINIMUM (Trough) to estimate envelope decay.
[minVal, idxMin] = min(Yo_fit);
tauMin = tau_fit(idxMin);

% Model at trough: |MinVal| approx c1 * exp(-c2 * tauMin^2)
% So: c2 = -log(|MinVal| / c1) / tauMin^2
if minVal < 0
    ratio = abs(minVal) / c_app(1);
    % Safety cap on ratio to prevent log errors
    ratio = max(ratio, 0.01); 
    c_app(2) = -log(ratio) / (tauMin^2);
else
    % If no negative dip found, fall back to a generic slow decay
    c_app(2) = 1 / (tau_fit(end)^2); 
end
c_app(2) = max(c_app(2), 1e-7);

% --- 9. Optimization with RELAXED BOUNDS ---
fun = @(c) (1./sqrt(Q_fit)) .* (Yo_fit - mCovF3(c, tau_fit));

% Lower bound: Amplitude > 0, Decay > 0, Freq can drop by 50%
lb = [c_app(1)*0.5;  0;  c_app(3)*0.5]; 

% Upper bound: Amp up to 150%, Decay up to 10x (allow fast drop), Freq up to 2x
ub = [c_app(1)*1.5;  c_app(2)*10;  c_app(3)*2.0];

options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'trust-region-reflective');
try
    coeffs = lsqnonlin(fun, c_app, lb, ub, options);
catch
    coeffs = c_app;
end
coeffs = coeffs(:);
end