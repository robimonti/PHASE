function [interp_values, interp_std, lonlatIN_AOI_NNI] = NaturalNeighborInterpolation(xyIN_AOI, values, std, ...
    xy_grid, centerline_data, figsDir, projDim, t_dateIN, lonlat_grid_AOI, utmZone, marker_size)

% NaturalNeighborInterpolation: Performs interpolation for each epoch and creates GIFs.
% Inputs:
%   xyIN_AOI:          - [N x 2] matrix of PS coordinates (UTM x, y).
%   values:            - [N x n_epochs] matrix of displacements for each PS and epoch.
%   std:               - [N x n_epochs] matrix of uncertainties for each PS and epoch.
%   xy_grid:           - [n_grid_points x 2] matrix of grid coordinates (UTM, 2D) or empty (1D).
%   centerline_data:   - struct with xy_centerline, s, projected_distances (1D) or empty (2D).
%   figsDir:           - directory to save GIF.
%   projDim:           - '1D' or '2D' for geometry type.
%   t_dateIN:          - [1 x n_epochs] datetime array of epochs.
%   lonlat_grid_AOI:   - [n_grid_points x 2] matrix of grid coordinates (lon, lat, 2D) or empty (1D).
%   utmZone:           - UTM zone string for coordinate conversion.
%   marker_size:       - marker size for the GIF plot
% Outputs:
%   interp_values: [n_grid_points x n_epochs] (2D) or [n_centerline_points x n_epochs] (1D).
%   interp_std: [n_grid_points x n_epochs] (2D) or [n_centerline_points x n_epochs] (1D).
%   lonlatIN_AOI_NNI


disp('--- Natural Neighbor Interpolation started... ---');


% Validate inputs
if isempty(xyIN_AOI) || isempty(values) || isempty(t_dateIN)
    error('xyIN_AOI, values, or t_dateIN is empty.');
end
n_epochs = size(values, 2);
if size(xyIN_AOI, 1) ~= size(values, 1)
    error('Number of PS points in xyIN_AOI (%d) does not match values (%d).', ...
        size(xyIN_AOI, 1), size(values, 1));
end

switch projDim
    case '2D'
        % Interpolate on precomputed grid
        grid_x = xy_grid(:,1);
        grid_y = xy_grid(:,2);
        interp_values = zeros(size(xy_grid, 1), n_epochs);   % displacement
        for t = 1:n_epochs
            F = scatteredInterpolant(xyIN_AOI(:,1), xyIN_AOI(:,2), values(:, t), 'natural', 'none');
            interp_values(:, t) = F(grid_x, grid_y);
        end
        interp_std = zeros(size(xy_grid, 1), n_epochs);   % uncertainty
        for t = 1:n_epochs
            F = scatteredInterpolant(xyIN_AOI(:,1), xyIN_AOI(:,2), std(:, t), 'natural', 'none');
            interp_std(:, t) = F(grid_x, grid_y);
        end

        % Create GIF with geobasemap
        fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
        v_min = min(round(prctile(values(:), 5), 0), -5); 
        v_max = max(round(prctile(values(:), 95), 0), 5);        
        [lat_ps, lon_ps] = utm2deg(xyIN_AOI(:,1), xyIN_AOI(:,2), repmat(utmZone, size(xyIN_AOI, 1), 1));
        h = waitbar(0, 'Plotting epochs...');
        for t = 1:n_epochs
            waitbar(t/n_epochs, h, sprintf('Plotting epoch %d/%d', t, n_epochs));
            geobasemap(gx, 'satellite')
            if t == 1
                pause(20)
            end
            hold(gx, 'on');
            geoscatter(gx,lonlat_grid_AOI(:,2), lonlat_grid_AOI(:,1), marker_size, interp_values(:, t), 'filled');
            geoscatter(gx,lat_ps, lon_ps, marker_size / 2, values(:, t), 'filled', 'MarkerEdgeColor', 'k');
            colormap(gx,jet);
            clim(gx,[v_min, v_max]);
            colorbar(gx);
            title(gx,sprintf('LOS displacement on %s (2D)', datestr(t_dateIN(t))), 'FontSize', 18);
            hold(gx, 'off');
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);
            if t == 1
                imwrite(imind, cm, fullfile(figsDir, 'NNI_displ2D.gif'), 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
            else
                imwrite(imind, cm, fullfile(figsDir, 'NNI_displ2D.gif'), 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
            end
            cla(gx);
        end
        close(fig);
        close(h);

        % Output
        lonlatIN_AOI_NNI = lonlat_grid_AOI;

    case '1D'
        % Validate centerline inputs
        if isempty(centerline_data) || isempty(centerline_data.xy_centerline) || isempty(centerline_data.s) || isempty(centerline_data.projected_distances)
            error('centerline_data or its fields are empty for 1D interpolation.');
        end

        % Interpolate along precomputed centerline
        grid_x = centerline_data.xy_centerline(:,1);
        grid_y = centerline_data.xy_centerline(:,2);
        s = centerline_data.s;
        projected_distances = centerline_data.projected_distances;
        interp_values = zeros(length(s), n_epochs);
        interp_std = zeros(length(s), n_epochs);

        % Interpolate along precomputed centerline using linear interpolation
        for t = 1:n_epochs
            % Remove invalid data
            valid_idx = ~isnan(projected_distances) & ~isnan(values(:, t));
            if sum(valid_idx) < 2
                warning('Epoch %d: Fewer than 2 valid points for 1D interpolation. Setting to NaN.', t);
                interp_values(:, t) = NaN;
                interp_std(:, t) = NaN;
                continue;
            end
            x = projected_distances(valid_idx); % PS points' projected distances
            v = values(valid_idx, t);           % displacement values for epoch t
            ss = std(valid_idx, t);             % displacement values for epoch t
            % Sort points and remove near-duplicates
            [x, sort_idx] = sort(x);
            v = v(sort_idx);
            ss = ss(sort_idx);
            unique_idx = [true; diff(x) > 1e-6]; % tolerance for near-duplicates
            x = x(unique_idx);
            v = v(unique_idx);
            ss = ss(unique_idx);
            if length(x) < 2
                warning('Epoch %d: Fewer than 2 unique projected distances (%d). Setting to NaN.', t, length(x));
                interp_values(:, t) = NaN;
                interp_std(:, t) = NaN;
                continue;
            end
            % Perform 1D linear interpolation directly
            try
                interp_values(:, t) = interp1(x, v, s, 'makima', NaN);
                interp_std(:, t) = interp1(x, ss, s, 'makima', NaN);
            catch e
                warning('Epoch %d: Interpolation failed (%s). Setting to NaN.', t, e.message);
                interp_values(:, t) = NaN;
                interp_std(:, t) = NaN;
            end
        end

        % 1D GIF with geobasemap
        fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
        [lat_centerline, lon_centerline] = utm2deg(grid_x, grid_y, repmat(utmZone, length(grid_x), 1));
        v_min = min(round(prctile(values(:), 5), 0), -5); 
        v_max = max(round(prctile(values(:), 95), 0), 5);
        [lat_ps, lon_ps] = utm2deg(xyIN_AOI(:,1), xyIN_AOI(:,2), repmat(utmZone, size(xyIN_AOI, 1), 1));
        h = waitbar(0, 'Plotting epochs...');
        for t = 1:n_epochs
            waitbar(t/n_epochs, h, sprintf('Plotting epoch %d/%d', t, n_epochs));
            geobasemap(gx, 'satellite')
            if t == 1
                pause(20)
            end
            hold(gx, 'on');
            geoscatter(gx,lat_centerline, lon_centerline, marker_size, interp_values(:, t), 'filled');
            geoscatter(gx,lat_ps, lon_ps, marker_size / 2, values(:, t), 'filled', 'MarkerEdgeColor', 'k');
            colormap(gx,jet);
            clim(gx,[v_min, v_max]);
            colorbar(gx);
            title(gx,sprintf('LOS displacement on %s (1D)', datestr(t_dateIN(t))), 'FontSize', 18);
            hold(gx, 'off');
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);
            if t == 1
                imwrite(imind, cm, fullfile(figsDir, 'NNI_displ1D.gif'), 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
            else
                imwrite(imind, cm, fullfile(figsDir, 'NNI_displ1D.gif'), 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
            end
            cla(gx);
        end
        close(fig);
        close(h);

        % Output
        lonlatIN_AOI_NNI = [lon_centerline, lat_centerline];
end


disp('=================================================');
disp('--- Natural Neighbor Interpolation started... ---');
disp('=================================================');