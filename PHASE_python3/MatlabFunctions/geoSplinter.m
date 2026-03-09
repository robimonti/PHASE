function [data, raster, par, sigmaPar] = geoSplinter (fileName, method);

% SINTAX:
%   [data, raster, par, sigmaPar] = geoSplinter (fileName, method);
%
% INPUT:
%   fileName = input file name (without extension)
%   method   = interpolation method:
%              - 'lin' = linear splines    (1D)
%              - 'cub' = cubic splines     (1D)
%              - 'bil' = bilinear splines  (2D)
%              - 'bic' = bicubic splines   (2D)
%
% OUTPUT:
%   data     = observation points, observed and estimated values
%   raster   = estimated values on a regular grid
%   par      = spline parameters
%   sigmaPar = standard deviations of the spline parameters
%
% DESCRIPTION:
%   Representation of the interpolated data.
%
% SEE ALSO   ...
%
% Mirko Reguzzoni
% OGS - Politecnico di Milano
% 22-01-2005

%-----------------------------------------------------------------------

if ((method == 'lin') | (method == 'cub'))
	dim = 1;
else % if ((method == 'bil') | (method == 'bic'))
	dim = 2;
end;

%-----------------------------------------------------------------------

data     = load ([fileName, '.out.txt']);
raster   = load ([fileName, '.grd.txt']);
par      = load ([fileName, '.par.txt']);
sigmaPar = load ([fileName, '.std.txt']);

%-----------------------------------------------------------------------

fid = fopen ([fileName, '.hdr.txt'], 'r');

str  = fscanf (fid, '%s %s', 2);
xMin = fscanf (fid, '%f', 1);

str  = fscanf (fid, '%s %s', 2);
xMax = fscanf (fid, '%f', 1);

if (dim == 2)

	str  = fscanf (fid, '%s %s', 2);
	yMin = fscanf (fid, '%f', 1);

	str  = fscanf (fid, '%s %s', 2);
	yMax = fscanf (fid, '%f', 1);

end;

str    = fscanf (fid, '%s %s', 2);
deltaX = fscanf (fid, '%f', 1);

xNum = round((xMax - xMin) / deltaX) + 1;
deltaX = (xMax - xMin) / (xNum - 1);

if (dim == 2)

	str    = fscanf (fid, '%s %s', 2);
	deltaY = fscanf (fid, '%f', 1);

	yNum = round((yMax - yMin) / deltaY) + 1;
	deltaY = (yMax - yMin) / (yNum - 1);

end;

% str    = fscanf (fid, '%s %s', 2);
% digits = fscanf (fid, '%f', 1);

fclose (fid);

%-----------------------------------------------------------------------

if (dim == 1)

	figure;
	subplot (2,1,1);
   plot (data(:,1), data(:,2), 'xb'); hold on;
   plot (data(:,1), data(:,3), 'xr'); hold off;
	title ('observed values (blue); estimated values (red)'); grid on;
	subplot (2,1,2);
	plot (data(:,1), data(:,4), 'xm-');
	title ('errors = observed values - estimated values'); grid on;

else % if (dim == 2)

	figure;
   stem3 (data(:,1), data(:,2), data(:,3), 'xb'); hold on;
   plot3 (data(:,1), data(:,2), data(:,4), 'xr'); hold off;
   xlabel('x-axis'); ylabel('y-axis');
   title ('observed values (blue); estimated values (red)'); grid on;

end;

%-----------------------------------------------------------------------

if (dim == 1)

  	figure;
   if (method == 'lin')
      plot (xMin:deltaX:xMax, raster, 'k');
   else % if (method == 'cub')
		xx = xMin : (xMax-xMin)/1000 : xMax;
		zz = interp1(xMin:deltaX:xMax, raster, xx, 'spline');
      plot (xx, zz, 'k');
   end;
   hold on;
   plot (data(:,1), data(:,2), 'xb');
   plot (data(:,1), data(:,3), 'xr');
   plot (xMin:deltaX:xMax, raster, 'xg'); hold off;
   title ('observed values (blue); estimated values (red); gridded values (green)');
   grid on;

else % if (dim == 2)

   figure;
   [xx yy] = meshgrid (xMin : (xMax-xMin)/100 : xMax, yMin : (yMax-yMin)/100 : yMax); 
   if (method == 'bil')
      zz = interp2(xMin:deltaX:xMax, yMin:deltaY:yMax, raster, xx, yy, 'linear');
   else % if (method == 'bic')
      zz = interp2(xMin:deltaX:xMax, yMin:deltaY:yMax, raster, xx, yy, 'spline');
   end;
   zzUD = flipud(zz);
   rasterUD = flipud(raster);
   surfl (xx, yy, zzUD); hold on;
   shading interp; colormap(pink);
   plot3 (data(:,1), data(:,2), data(:,3), 'xb');
   plot3 (data(:,1), data(:,2), data(:,4), 'xr');
   stem3 (xMin:deltaX:xMax, yMin:deltaY:yMax, rasterUD, 'xg'); hold off;
   xlabel('x-axis'); ylabel('y-axis');
   title ('observed values (blue); estimated values (red); gridded values (green)');
   grid on;

   figure;
   contour (xx, yy, zzUD); colorbar; grid on;
   xlabel('x-axis'); ylabel('y-axis');
   title ('interpolated surface (contour)');

   figure;
	mesh (xMin:deltaX:xMax, yMin:deltaY:yMax, rasterUD); colorbar;
   xlabel('x-axis'); ylabel('y-axis');
   title ('gridded values (mesh)');

   figure;
	imagesc (raster); colorbar;
	axis  ('square');
   title ('gridded values (image)');

end;

%-----------------------------------------------------------------------

if (dim == 1)

	figure;
	errorbar (xMin:deltaX:xMax, par, sigmaPar, 'xb'); grid on;
   title ('spline parameters and their standard deviations');

else % if (dim == 2)

   figure;
	imagesc (par); colorbar;
	axis  ('square');
   title ('spline parameters (image)');

   figure;
	imagesc (sigmaPar); colorbar;
	axis  ('square');
   title ('parameter std. dev. (image)');

end;

%-----------------------------------------------------------------------