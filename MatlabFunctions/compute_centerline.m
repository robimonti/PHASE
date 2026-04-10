function [centerline_data, xy_grid, lonlat_grid_AOI] = compute_centerline(xyAOI, cline_resolution, xyIN_AOI)

% compute_centerline computes the centerline of a polygon for 1D geospatial analysis.
%
% Inputs:
%   xyAOI            - Nx2 matrix of polygon vertices in UTM coordinates [x, y].
%   cline_resolution - Spatial resolution for centerline interpolation [meters].
%   xyIN_AOI         - Mx2 matrix of PS points inside AOI in UTM coordinates [x, y].
%
% Outputs:
%   centerline_data  - Struct containing:
%                      - xy_centerline: Interpolated centerline coordinates [x, y].
%                      - s: Cumulative distances along interpolated centerline.
%                      - projected_distances: Distances of PS points projected onto centerline.
%   xy_grid          - Empty array (for 2D case compatibility).
%   lonlat_grid_AOI  - Empty array (for 2D case compatibility).
%
% Notes:
%   - Requires Image Processing Toolbox for bwmorph.
%   - Assumes xyAOI is a closed polygon; closes it if not.
%   - Pixel size is set to cline_resolution/5 for smooth skeletonization.

    % Validate inputs
    if size(xyAOI, 1) < 3
        error('xyAOI must have at least 3 points to define a polygon.');
    end
    if cline_resolution <= 0
        error('cline_resolution must be positive.');
    end

    % Ensure xyAOI is strictly unique vertices (remove closing duplicate if it exists)
    if isequal(xyAOI(1,:), xyAOI(end,:))
        xyAOI_unique = xyAOI(1:end-1, :);
    else
        xyAOI_unique = xyAOI;
    end

    num_vertices = size(xyAOI_unique, 1);

    % =====================================================================
    % ENGINE 1: EXACT ANALYTICAL SOLUTION FOR 4-VERTEX POLYGONS (BRIDGES)
    % =====================================================================
    if num_vertices == 4
        % Edges of the quadrilateral
        L12 = norm(xyAOI_unique(2,:) - xyAOI_unique(1,:));
        L23 = norm(xyAOI_unique(3,:) - xyAOI_unique(2,:));
        L34 = norm(xyAOI_unique(4,:) - xyAOI_unique(3,:));
        L41 = norm(xyAOI_unique(1,:) - xyAOI_unique(4,:));
        
        % Identify which pair of opposing edges are the "short" ends
        avg_pair_A = (L12 + L34) / 2;
        avg_pair_B = (L23 + L41) / 2;
        
        if avg_pair_A < avg_pair_B
            pt1 = (xyAOI_unique(1,:) + xyAOI_unique(2,:)) / 2;
            pt2 = (xyAOI_unique(3,:) + xyAOI_unique(4,:)) / 2;
        else
            pt1 = (xyAOI_unique(2,:) + xyAOI_unique(3,:)) / 2;
            pt2 = (xyAOI_unique(4,:) + xyAOI_unique(1,:)) / 2;
        end
        
        % Force start_pt to be the southernmost point (minimum Y)
        if pt1(2) < pt2(2)
            start_pt = pt1;
            end_pt = pt2;
        else
            start_pt = pt2;
            end_pt = pt1;
        end
        
        % Generate interpolated centerline
        L_total = norm(end_pt - start_pt);
        s_sampled = (0:cline_resolution:L_total)';
        if s_sampled(end) ~= L_total
            s_sampled = [s_sampled; L_total];
        end
        
        dir_vec = (end_pt - start_pt) / L_total;
        xy_centerline = [start_pt(1) + s_sampled * dir_vec(1), start_pt(2) + s_sampled * dir_vec(2)];
        s = s_sampled;

    % =====================================================================
    % ENGINE 2: UPGRADED SKELETONIZATION FOR COMPLEX CURVES
    % =====================================================================
    else
        % Close the polygon for masking
        xyAOI_closed = [xyAOI_unique; xyAOI_unique(1,:)];
        
        x_min = min(xyAOI_closed(:,1)); x_max = max(xyAOI_closed(:,1));
        y_min = min(xyAOI_closed(:,2)); y_max = max(xyAOI_closed(:,2));

        % Dynamic high-resolution pixel size (Max 1 meter to preserve thin bridges)
        pixel_size = min(cline_resolution / 4, 1.0); 
        nx = ceil((x_max - x_min) / pixel_size);
        ny = ceil((y_max - y_min) / pixel_size);

        % Create binary image
        [X, Y] = meshgrid(linspace(x_min, x_max, nx), linspace(y_min, y_max, ny));
        binary_image = inpolygon(X, Y, xyAOI_closed(:,1), xyAOI_closed(:,2));

        % Skeletonize and clean tiny branches (spurs)
        skeleton = bwmorph(binary_image, 'skel', Inf);
        skeleton = bwmorph(skeleton, 'spur', 5);

        % Find endpoints to start tracing
        eps = bwmorph(skeleton, 'endpoints');
        [ep_y, ep_x] = find(eps);
        [skel_y, skel_x] = find(skeleton);

        if isempty(skel_x)
            error('Skeletonization failed. AOI might be too narrow or degenerate.');
        end

        % Force tracing to start from the southernmost endpoint
        start_idx = 1;
        if ~isempty(ep_y)
            % Convert pixel row indices to UTM Y coordinates
            ep_utm_y = y_min + (ep_y - 1) * pixel_size;
            % Find the index of the endpoint with the minimum Y
            [~, south_idx] = min(ep_utm_y); 
            % Set the starting node for the graph tracer
            start_idx = find(skel_x == ep_x(south_idx) & skel_y == ep_y(south_idx), 1);
        end

        % Trace the skeleton sequentially
        num_skel = length(skel_x);
        ordered_idx = zeros(num_skel, 1);
        visited = false(num_skel, 1);
        curr_idx = start_idx;

        for k = 1:num_skel
            visited(curr_idx) = true;
            ordered_idx(k) = curr_idx;
            unvisited = find(~visited);
            if isempty(unvisited), break; end
            
            % Find nearest unvisited neighbor
            dist_sq = (skel_x(unvisited) - skel_x(curr_idx)).^2 + (skel_y(unvisited) - skel_y(curr_idx)).^2;
            [min_dist_sq, closest_unvisited] = min(dist_sq);
            
            if min_dist_sq > 9 % If gap is larger than 3 pixels, curve is broken
                break;
            end
            curr_idx = unvisited(closest_unvisited);
        end

        % Clean unused pre-allocated zeros and convert to UTM
        ordered_idx = ordered_idx(ordered_idx > 0);
        skel_utm_x = x_min + (skel_x(ordered_idx) - 1) * pixel_size;
        skel_utm_y = y_min + (skel_y(ordered_idx) - 1) * pixel_size;
        xy_skeleton_ordered = [skel_utm_x, skel_utm_y];

        % Calculate cumulative distances for interpolation
        dist_steps = sqrt(diff(xy_skeleton_ordered(:,1)).^2 + diff(xy_skeleton_ordered(:,2)).^2);
        
        % Remove duplicate consecutive points for stable interpolation
        valid_idx = [true; dist_steps > 0];
        xy_skeleton_ordered = xy_skeleton_ordered(valid_idx, :);
        cumdist = [0; cumsum(dist_steps(dist_steps > 0))];
        
        s_total = cumdist(end);
        s = (0:cline_resolution:s_total)';
        if s(end) ~= s_total, s = [s; s_total]; end
        
        centerline_x = interp1(cumdist, xy_skeleton_ordered(:,1), s, 'linear');
        centerline_y = interp1(cumdist, xy_skeleton_ordered(:,2), s, 'linear');
        xy_centerline = [centerline_x, centerline_y];
    end

    % =====================================================================
    % STEP 3: MATHEMATICAL PROJECTION OF PS POINTS ONTO CENTERLINE
    % =====================================================================
    num_ps = size(xyIN_AOI, 1);
    projected_distances = zeros(num_ps, 1);
    
    for i = 1:num_ps
        P = xyIN_AOI(i, :);
        min_dist = Inf;
        proj_s = 0;
        
        for j = 1:size(xy_centerline, 1)-1
            A = xy_centerline(j, :);
            B = xy_centerline(j+1, :);
            AB = B - A;
            AP = P - A;
            
            norm_AB_sq = dot(AB, AB);
            if norm_AB_sq == 0, continue; end
            
            t = dot(AP, AB) / norm_AB_sq;
            t = max(0, min(1, t)); % Clamp to segment
            
            proj_point = A + t * AB;
            dist_to_segment = norm(P - proj_point);
            
            if dist_to_segment < min_dist
                min_dist = dist_to_segment;
                seg_s_start = s(j);
                seg_length = norm(AB);
                proj_s = seg_s_start + t * seg_length;
            end
        end
        projected_distances(i) = proj_s;
    end

    % Prepare Outputs
    centerline_data = struct('xy_centerline', xy_centerline, 's', s, 'projected_distances', projected_distances);
    xy_grid = [];
    lonlat_grid_AOI = [];
    
    disp('Robust centerline projection complete.');
end