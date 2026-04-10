function [accum_matrix, sum_sq, n_pts] = f2DEmpCovEst_accum(S, Coords, bins)

% f2DEmpCovEst_accum Estimates Empirical Covariance components for 2D spatial data.
%
% Inputs:
%   S      : Signal vector (N x 1)
%   Coords : Spatial coordinates matrix (N x 2) -> [x, y]
%   bins   : Vector defining the edges of the distance bins
%
% Outputs:
%   accum_matrix : [Sum_Products, Count, Sum_Distances] per bin (Lags > 0)
%   sum_sq       : Sum of S.^2 (for calculating global Variance at Lag 0)
%   n_pts        : Number of points (count for Lag 0)

S = S(:);
nS = length(S);

% --- 1. Variance Calculation (Lag 0) ---
sum_sq = sum(S.^2);
n_pts  = nS;

% Initialize accumulator
n_bins = length(bins) - 1;
accum_matrix = zeros(n_bins, 3);

% If fewer than 2 points, no pairs can be formed
if nS < 2, return; end

% --- 2. Covariance Calculation (Lag > 0) ---
% Generate indices for the upper triangle (unique pairs i, j where i < j)
[I, J] = find(triu(true(nS), 1));

if isempty(I), return; end

% Calculate Products
prod_vals = S(I) .* S(J);

% Calculate 2D Euclidean Distances
% d = sqrt( (x1-x2)^2 + (y1-y2)^2 )
d_x = Coords(I, 1) - Coords(J, 1);
d_y = Coords(I, 2) - Coords(J, 2);
dist_vals = sqrt(d_x.^2 + d_y.^2);

% --- BINNING LOGIC ---
% Discretize distances into defined bins
bin_idx = discretize(dist_vals, bins);

% Filter out pairs that fall outside the defined bins (NaNs)
valid = ~isnan(bin_idx);

if ~any(valid), return; end

b = bin_idx(valid);
p = prod_vals(valid);
d = dist_vals(valid);

% Accumulate results
% Column 1: Sum of Products (Covariance numerator)
% Column 2: Count of pairs
% Column 3: Sum of Distances (Mean lag calculation)
accum_matrix = [accumarray(b, p, [n_bins, 1]), ...
                accumarray(b, 1, [n_bins, 1]), ...
                accumarray(b, d, [n_bins, 1])];
end