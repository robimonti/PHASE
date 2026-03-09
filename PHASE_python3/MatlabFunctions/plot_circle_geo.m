function [lonlat_circle] = plot_circle_geo(centre_lon, centre_lat, radius)

% plot_circle_geo Plots a circle on a geographic plot given center coordinates and radius.
%
% Input:
%   centre_lon: Longitude of the circle's center (degrees).
%   centre_lat: Latitude of the circle's center (degrees).
%   radius: Radius of the circle in meters.
% Output:
%   lonlat_circle: matrix containing longitude and latitude of circle points

% Earth's radius (mean radius) in meters
R_earth = 6371000;

% Convert radius to angular distance (radians)
angular_radius = radius / R_earth;

% Number of points to define the circle
num_points = 100;

% Generate points around the circle
theta = linspace(0, 2*pi, num_points);
lat_rad = deg2rad(centre_lat);
lon_rad = deg2rad(centre_lon);

lat_circle = asin(sin(lat_rad) * cos(angular_radius) + ...
                  cos(lat_rad) * sin(angular_radius) * cos(theta));
lon_circle = lon_rad + atan2(sin(theta) * sin(angular_radius) * cos(lat_rad), ...
                             cos(angular_radius) - sin(lat_rad) * sin(lat_circle));

% Convert back to degrees
lat_circle = rad2deg(lat_circle);
lon_circle = rad2deg(lon_circle);

lonlat_circle = [lon_circle', lat_circle'];

end