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

% Ensure xyAOI is closed
if ~isequal(xyAOI(1,:), xyAOI(end,:))
    xyAOI = [xyAOI; xyAOI(1,:)];
end

% - 1: Convert polygon to a binary image for skeletonization
x_min = min(xyAOI(:,1));
x_max = max(xyAOI(:,1));
y_min = min(xyAOI(:,2));
y_max = max(xyAOI(:,2));

% Define image resolution (meters per pixel)
pixel_size = 10; % [m]
nx = ceil((x_max - x_min) / pixel_size);
ny = ceil((y_max - y_min) / pixel_size);

% Create binary image of the polygon
[X, Y] = meshgrid(linspace(x_min, x_max, nx), linspace(y_min, y_max, ny));
points = [X(:), Y(:)];
in_poly = inpolygon(points(:,1), points(:,2), xyAOI(:,1), xyAOI(:,2));
binary_image = reshape(in_poly, ny, nx);

% - 2: Apply morphological skeletonization
skeleton = bwmorph(binary_image, 'skel', Inf);

% - 3: Convert skeleton back to vector coordinates
[skel_y, skel_x] = find(skeleton);
if isempty(skel_x)
    error('Skeletonization failed to produce a centerline.');
end
skel_x = x_min + (skel_x - 1) * pixel_size;
skel_y = y_min + (skel_y - 1) * pixel_size;
xy_skeleton = [skel_x, skel_y];

% - 4: Simplify skeleton to a single path, starting from southernmost point
% Find southernmost point (minimum y-coordinate)
[~, start_idx] = min(xy_skeleton(:,2));
current_point = xy_skeleton(start_idx, :);
xy_skeleton(start_idx, :) = [];
xy_centerline = current_point;

while ~isempty(xy_skeleton)
    distances = sqrt((xy_skeleton(:,1) - current_point(1)).^2 + (xy_skeleton(:,2) - current_point(2)).^2);
    [min_dist, idx] = min(distances);
    if min_dist > pixel_size * 2 % Break if points are too far apart
        break;
    end
    current_point = xy_skeleton(idx, :);
    xy_centerline = [xy_centerline; current_point];
    xy_skeleton(idx, :) = [];
end

% Ensure centerline is oriented from south to north
[~, sort_idx] = sort(xy_centerline(:,2), 'ascend');
xy_centerline = xy_centerline(sort_idx, :);

% - 5: Interpolate centerline at cline_resolution intervals
dist = sqrt(diff(xy_centerline(:,1)).^2 + diff(xy_centerline(:,2)).^2);
cumdist = [0; cumsum(dist)];
s_total = cumdist(end);
if s_total == 0
    error('Centerline has zero length. Check skeletonization output.');
end
s = 0:cline_resolution:s_total;
centerline_x = interp1(cumdist, xy_centerline(:,1), s, 'linear');
centerline_y = interp1(cumdist, xy_centerline(:,2), s, 'linear');
xy_centerline = [centerline_x', centerline_y'];

% - 6: Project PS points onto centerline segments
projected_distances = zeros(size(xyIN_AOI, 1), 1);
for i = 1:size(xyIN_AOI, 1)
    P = xyIN_AOI(i, :); % PS point
    min_dist = Inf;
    proj_s = 0;
    
    % Check each segment of the interpolated centerline
    for j = 1:size(xy_centerline, 1)-1
        A = xy_centerline(j, :); % Segment start
        B = xy_centerline(j+1, :); % Segment end
        AB = B - A; % Segment vector
        AP = P - A; % Vector from segment start to PS point
        
        % Compute projection parameter t (0 <= t <= 1 for within segment)
        t = dot(AP, AB) / dot(AB, AB);
        t = max(0, min(1, t)); % Clamp to segment endpoints
        
        % Compute projection point
        proj_point = A + t * AB;
        
        % Compute distance from PS point to projection point
        dist = norm(P - proj_point);
        
        % Update if this is the closest projection
        if dist < min_dist
            min_dist = dist;
            % Interpolate cumulative distance for projection point
            seg_length = norm(AB);
            seg_s_start = s(j);
            seg_s_end = s(j+1);
            proj_s = seg_s_start + t * (seg_s_end - seg_s_start);
        end
    end
    projected_distances(i) = proj_s;
end

% - 7: Store outputs
centerline_data = struct('xy_centerline', xy_centerline, 's', s', 'projected_distances', projected_distances);
xy_grid = [];
lonlat_grid_AOI = [];

end