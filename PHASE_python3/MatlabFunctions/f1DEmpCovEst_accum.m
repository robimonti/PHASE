function [accum_matrix, sum_sq, n_pts] = f1DEmpCovEst_accum(S, t, bins, is_space)

% f1DEmpCovEst_accum 
%
% Outputs:
%   accum_matrix : For Lags > 0
%   sum_sq       : Sum of S.^2 (for calculating global Variance at Lag 0)
%   n_pts        : Number of points (count for Lag 0)

S = S(:); 
t = t(:);
nS = length(S);

% --- 1. Variance Calculation (Lag 0) ---
% This is strictly the diagonal elements
sum_sq = sum(S.^2);
n_pts  = nS;

% Output initialization for Covariance
if is_space
    n_bins = length(bins) - 1;
    accum_matrix = zeros(n_bins, 3);
else
    dt = bins;
    accum_matrix = []; % Dynamic
end

if nS < 2, return; end

% --- 2. Covariance Calculation (Lag > 0) ---
[I, J] = find(triu(true(nS), 1)); 

if isempty(I), return; end

prod_vals = S(I) .* S(J);
dist_vals = abs(t(I) - t(J));

% --- BINNING LOGIC ---
if is_space
    bin_idx = discretize(dist_vals, bins);
    valid = ~isnan(bin_idx);
    if ~any(valid), return; end
    
    b = bin_idx(valid);
    p = prod_vals(valid);
    d = dist_vals(valid);
    
    accum_matrix = [accumarray(b, p, [n_bins, 1]), ...
                    accumarray(b, 1, [n_bins, 1]), ...
                    accumarray(b, d, [n_bins, 1])];
else
    dt = bins;
    bin_idx = round(dist_vals / dt); % Index 1 = Lag 1*dt (since dist > 0)
    
    % Safety: ensure bin_idx >= 1
    bin_idx(bin_idx < 1) = 1; 
    
    max_lag_found = max(bin_idx);
    accum_matrix = [accumarray(bin_idx, prod_vals, [max_lag_found, 1]), ...
                    accumarray(bin_idx, 1, [max_lag_found, 1]), ...
                    accumarray(bin_idx, dist_vals, [max_lag_found, 1])];
end
end