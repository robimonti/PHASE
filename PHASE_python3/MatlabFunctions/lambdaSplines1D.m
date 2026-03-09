function lambda = lambdaSplines1D(observations, sigma2v, tauGrid, method)
    % Computes the regularization parameter for 1D splines interpolation
    %
    % Inputs:
    % - observations: Nx2 matrix containing (x, z) data points
    % - sigma2v: variance of the noise
    % - tauGrid: step of the splines grid
    % - method: 1 for linear or 2 for cubic
    %
    % Output:
    % - lambda: regularization parameter
    %
    % (c) Roberto Monti
    % Politecnico di Milano
    % Last update: Jan. 2025

    if nargin < 4
        error('Please specify the interpolation method: "1 - linear" or "2 - cubic".');
    end

    % Ensure the input is in the correct format
    if size(observations, 2) ~= 2
        error('Points must be an Nx2 matrix with columns (x, z).');
    end

    % Extract coordinates
    x = observations(:, 1);
    z = observations(:, 2);

    % Number of points
    numPoints = size(observations, 1);

    % Initialize derivatives
    derivatives = [];

    % Compute derivatives for each point
    for i = 1:numPoints
        % Find the closest point
        distances = abs(x - x(i));
        distances(i) = inf; % Exclude self
        [~, closestIdx] = min(distances);

        % Compute first derivative
        dist = distances(closestIdx);
        if dist == 0
            continue; % Skip if distance is zero
        end
        dz = z(closestIdx) - z(i);
        derivative = dz / dist;

        % Store derivative
        derivatives = [derivatives; derivative];
    end

    if method == 2   % cubic
        % Compute second derivatives
        secondDerivatives = [];
        for i = 1:numel(derivatives) - 1
            dDeriv = derivatives(i + 1) - derivatives(i);
            secondDerivative = dDeriv / (x(i + 1) - x(i));
            secondDerivatives = [secondDerivatives; secondDerivative];
        end

        % Update derivatives for sigma2n calculation
        derivatives = secondDerivatives;
    end

    % Compute sigma^2_n
    sigma2n = sum(derivatives.^2) / numel(derivatives);

    % Compute lambda
    if method == 1
        lambda = sigma2v / (sigma2n * tauGrid^2);
    else
        lambda = sigma2v / (sigma2n * tauGrid^4);
    end
end