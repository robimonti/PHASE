function lambda = lambdaSplines2D(observations, sigma2v, deltaGrid, method)
    % Computes the regularization parameter for splines interpolation
    % 
    % Inputs:
    % - observations: Nx3 matrix containing (x, y, z) data points
    % - sigma2v: variance of the noise
    % - deltaGrid: step of the splines grid
    % - method: 1 for bilinear or 2 for bicubic
    %
    % Output:
    % - lambda: regularization parameter
    %
    %
    % (c) Roberto Monti
    % Politecnico di Milano
    % Last update: Dec. 2024

    
    if nargin < 4
        error('Please specify the interpolation method: "1 - bilinear" or "2 - bicubic".');
    end

    % Ensure the input is in the correct format
    if size(observations, 2) ~= 3
        error('Points must be an Nx3 matrix with columns (x, y, z).');
    end

    % Extract coordinates
    x = observations(:, 1);
    y = observations(:, 2);
    z = observations(:, 3);
    
    % Number of points
    numPoints = size(observations, 1);
    
    % Initialize derivatives and corresponding midpoints
    derivatives = [];
    midpoints = [];
    dist1 = [];
    
    % Loop over each point to compute the derivatives
    for i = 1:numPoints
        % Find the closest point
        distances = sqrt((x - x(i)).^2 + (y - y(i)).^2);
        distances(i) = inf; % Exclude self
        [~, closestIdx] = min(distances);
        
        % Compute first derivative
        dist = distances(closestIdx);
        if dist == 0
            continue; % Skip if distance is zero
        end
        dz = z(closestIdx) - z(i);
        derivative = dz / dist;
        
        % Compute the midpoint
        midpoint = [(x(i) + x(closestIdx)) / 2, ...
                    (y(i) + y(closestIdx)) / 2];
        
        % Store derivative and midpoint
        derivatives = [derivatives; derivative];
        midpoints = [midpoints; midpoint];
        dist1 = [dist1;dist];
    end

    if method == 2   % bicubic
        % Remove duplicate midpoints (keep one derivative for each unique midpoint)
        [uniqueMidpoints, uniqueIdx] = unique(midpoints, 'rows');
        derivatives = derivatives(uniqueIdx);

        % Recompute derivatives for the new coordinates
        newDerivatives = [];
        newMidpoints = [];
        numMidpoints = size(uniqueMidpoints, 1);
        for i = 1:numMidpoints
            distances = sqrt(sum((uniqueMidpoints - uniqueMidpoints(i, :)).^2, 2));
            distances(i) = inf;
            [~, closestIdx] = min(distances);
            
            dist = distances(closestIdx);
            if dist == 0
                continue; % Skip if distance is zero
            end
            dDeriv = derivatives(closestIdx) - derivatives(i);
            newDerivative = dDeriv / dist;
            midpoint = (uniqueMidpoints(i, :) + uniqueMidpoints(closestIdx, :)) / 2;

            % Store new derivative and midpoint
            newDerivatives = [newDerivatives; newDerivative];
            newMidpoints = [newMidpoints; midpoint];
        end

        % Update derivatives for sigma2n calculation
        derivatives = newDerivatives;
    end

    % Compute sigma^2_n
    sigma2n = sum(derivatives.^2) / numel(derivatives);

    % Compute lambda
    if method == 1
        lambda = sigma2v / (sigma2n * deltaGrid^2);
    else
        lambda = sigma2v / (sigma2n * deltaGrid^4);
    end

end