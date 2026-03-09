function [sigma, eCov, sigmaGrid, eCovF, Cecf] = f1DEmpCovEst_fast(S, t, dsCov)

% f1DEmpCovEst_fast Optimized 1D Empirical Covariance Estimation

nS = length(S);

% Get all pairs only once (upper triangular, excluding diag)
[i1, i2] = find(triu(true(nS), 1));

% Distances and products
sigma = abs(t(i1) - t(i2));
eCov  = S(i1) .* S(i2);

% Bin edges
sigmaMax = max(sigma);
sigmaLim = (0:dsCov:sigmaMax)';
nBins = length(sigmaLim);

% Bin assignment
[~,~,bin] = histcounts(sigma, sigmaLim);

% Prepare outputs
sigmaGrid = zeros(nBins-1,1);
eCovF     = zeros(nBins-1,1);
Cecf      = zeros(nBins-1,1);

for k = 1:nBins-1
    idx = (bin == k);
    if any(idx)
        sigmaGrid(k) = mean(sigma(idx));
        eCovF(k)     = mean(eCov(idx));
        Cecf(k)      = 1 / sum(idx);
    else
        sigmaGrid(k) = NaN;
        eCovF(k)     = NaN;
        Cecf(k)      = NaN;
    end
end

% Handle sigma=0 separately
sigma = [0; sigma(:)];
eCov  = [mean(S.^2); eCov(:)];
sigmaGrid = [0; sigmaGrid];
eCovF     = [mean(S.^2); eCovF];
Cecf      = [(1/nS)^2; Cecf];

% Remove NaNs (in case some bins are empty)
idx = isnan(sigmaGrid) | isnan(eCovF) | isnan(Cecf);
sigmaGrid(idx) = [];
eCovF(idx) = [];
Cecf(idx) = [];

end