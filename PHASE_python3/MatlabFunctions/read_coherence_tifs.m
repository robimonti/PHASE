function coh_values = read_coherence_tifs(xyIN_AOI, t, coherence_dir, master_date, utmZone)

% read_coherence_tifs Reads SAR coherence for PS
%
% Inputs:
%   xyIN_AOI: [N x 2] matrix of PS UTM coordinates [x, y]
%   t: [M x 1] vector of epoch timestamps (days since master)
%   coherence_dir: Path to folder with coherence .tif files
%   master_date: Datetime of master epoch (e.g., datetime('2019-07-03'))
%   utmZone: [N x 4] matrix of UTM zones (e.g., '33 N')
% Output:
%   coh_values: [N x M] matrix of coherence values per PS and epoch

% Initialize
N_ps = size(xyIN_AOI, 1);
M = length(t);
coh_values = nan(N_ps, M);

% List .tif files
tif_files = dir(fullfile(coherence_dir, '*_IW*_coh_TC.tif'));
if isempty(tif_files)
    error('No coherence .tif files found in %s', coherence_dir);
end

% Convert t to datetimes
t_dates = master_date + days(t);

% Verify UTM zone consistency
unique_zones = unique(utmZone, 'rows');
if size(unique_zones, 1) > 1
    warning('Multiple UTM zones detected: %s. Using first zone: %s', ...
            strjoin(unique_zones, ', '), unique_zones(1, :));
end
target_zone = unique_zones(1, :);

% Process each .tif file
for k = 1:length(tif_files)
    % Extract master and slave dates
    fname = tif_files(k).name;
    date_tokens = regexp(fname, '(\d{8})_(\d{8})_IW\d+_coh_TC\.tif', 'tokens');
    if isempty(date_tokens)
        warning('Filename %s does not match expected format.', fname);
        continue;
    end
    dates = date_tokens{1};
    slave_date = datetime(dates{2}, 'InputFormat', 'yyyyMMdd');

    % Match to epoch
    [~, idx_t] = min(abs(t_dates - slave_date));
    if abs(t_dates(idx_t) - slave_date) > days(1)
        continue; % Skip if no close match
    end

    % Read GeoTIFF
    tif_file = fullfile(coherence_dir, fname);
    [coh_img, R] = geotiffread(tif_file);
    info = geotiffinfo(tif_file);

    % Verify .tif UTM zone
    if isfield(info.GeoTIFFTags.GeoKeyDirectoryTag, 'ProjectedCSTypeGeoKey')
        tif_zone_num = info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey - 32600;
        if info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey >= 32601 && ...
           info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey <= 32660
            tif_zone_letter = 'N';
        else
            tif_zone_letter = 'S';
        end
        tif_zone = sprintf('%d %c', tif_zone_num, tif_zone_letter);
    else
        error('ProjectedCSTypeGeoKey missing in %s.', fname);
    end
    if ~strcmp(tif_zone, target_zone)
        warning('UTM zone mismatch: .tif (%s) vs PS (%s) for %s', tif_zone, target_zone, fname);
        continue;
    end

    % Extract coherence at PS locations
    for i = 1:N_ps
        % Convert UTM x,y to pixel coordinates
        [col, row] = worldToIntrinsic(R, xyIN_AOI(i, 1), xyIN_AOI(i, 2));
        row = round(row);
        col = round(col);

        % Check bounds and extract
        if row >= 1 && row <= size(coh_img, 1) && col >= 1 && col <= size(coh_img, 2)
            coh = coh_img(row, col);
            if coh >= 0 && coh <= 1
                coh_values(i, idx_t) = coh;
            end
        end
    end
end

% Check pixel spacing
if ~isempty(tif_files)
    info = geotiffinfo(fullfile(coherence_dir, tif_files(1).name));
    % Hardcode pixel spacing from QGIS
    pixel_spacing = [11.02413963681841658, 14.12940703933526443];
    % Attempt to parse metadata
    try
        if isfield(info, 'GeoTIFFTags') && isfield(info.GeoTIFFTags, 'ModelPixelScaleTag')
            pixel_spacing = [info.GeoTIFFTags.ModelPixelScaleTag(1), info.GeoTIFFTags.ModelPixelScaleTag(2)];
        elseif isfield(info, 'GeoTIFFTags') && isfield(info.GeoTIFFTags, 'ModelTransformationTag')
            pixel_spacing = [info.GeoTIFFTags.ModelTransformationTag(1, 1), abs(info.GeoTIFFTags.ModelTransformationTag(2, 2))];
        else
            warning('Pixel spacing metadata missing in %s. Using expected values: %.2f x %.2f m.', ...
                    tif_files(1).name, pixel_spacing(1), pixel_spacing(2));
        end
    catch
        warning('Error parsing pixel spacing metadata in %s. Using expected values: %.2f x %.2f m.', ...
                tif_files(1).name, pixel_spacing(1), pixel_spacing(2));
    end
    fprintf('Coherence .tif pixel spacing: %.2f x %.2f m\n', pixel_spacing);
end
end