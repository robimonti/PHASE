function [sigma, eCov, sigmaGrid, eCovF, Cecf, h] = f1DEmpCovEst(S, t, dsCov, fig)
% function to compute Empirical Covariance function
%
% SYNTAX: [sigma, eCov, sigmaGrid, eCovF, Cecf, h1, h2] = fCartCovEmpEst(S, x, y, dCov)
%
% INPUT:
% - S        --> Signal for which the covariance will be estimated
% - t        --> coordinate vector
% - dsCov    --> step for averaging the point cloud of covariances
% - fig      --> 0 --> no figure,
%                1 --> only averaged covariance,
%                2 --> averaged + sparse points
%
% OUTPUT:
% - sigmaGrid --> sample distance of empirical covariance function
% - eCovF     --> empirical covariance function sampled at spherical distance sigmaGrid
% - Cecf      --> diagonal of the cofactor matrix of empirical covariance function
% - h         --> handle of the figures (size depending on fig input)

h = -1;  % initialize the handle to the figures

% number of signal observations
nS = length(S);

% combination of coordinates of all the points
T = repmat(t, 1, nS);

% compute all the distances
sigma = abs(T - T');
clear T

% compute the covariance for all the possible couples of points
Sall = repmat(S, 1, nS);
eCov = Sall .* Sall';
clear Sall

% define the limits of the classes in term of psi values
sigmaMax = max(sigma(:));
sigmaLim = (0:dsCov:sigmaMax)';

% intitalize output variables
eCovF = zeros(length(sigmaLim), 1);      % empirical covariance
Cecf  = zeros(length(sigmaLim), 1);      % cofactor matrix (proportional to the inverse of number of data of each class)
sigmaGrid = zeros(length(sigmaLim), 1);  % grid of distance at which covariance will be computed

% distance equal to 0
% eCovF(1)  = mean(eCov(sigma(:)==0)) * sum(eCov(sigma(:)==0)) / (sum(eCov(sigma(:)==0))-1);
eCovF(1)  = mean(eCov(sigma(:)==0));
Cecf(1) = (1./numel(eCov(sigma==0))).^2;

% starting from 2, since the class 1 is the zero distance
for i = 2:length(sigmaLim) 
    idx = (sigma > sigmaLim(i-1)) & (sigma <= sigmaLim(i));  % points inside the current class
    % eCovF(i)     = sum(eCov(idx))./(sum(idx(:)) - i);       % empirical covariance of the class
    eCovF(i)     = mean(eCov(idx));
    sigmaGrid(i) = mean(sigma(idx));      % mean sigma of the class (to be used for interpolating the model)
    % sigmaGrid(i) = (sigmaLim(i-1) + sigmaLim(i))/2;
    Cecf(i)      = 1./sum(idx(:));        % the cofactor is proportional to the inverse of the number of data
end
idx = find(isnan(eCovF) | isnan(sigmaGrid) | isnan(Cecf));
sigmaGrid(idx) = [];
eCovF(idx)     = [];
Cecf(idx)      = [];

if fig >= 1
    h(1) = figure; plot(sigmaGrid,eCovF,'.');
    xlim([0 max(sigmaGrid)/2]);
    xlabel('distance');
    title('Empirical covariance function');
end
if fig >= 2
    h(2) = figure;
    plot(sigma(:), eCov(:), '.')
    xlabel('distance');
    title('Covariance point cloud');
end
sigma = sigma(:);
eCov = eCov(:);