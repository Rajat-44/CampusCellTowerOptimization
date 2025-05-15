clc; clear; close all;

grid_width = 1400;
grid_height = 540;

lat_min = 29.861973;  lat_max = 29.864156;
lon_min = 77.895126;  lon_max = 77.901626;

lat_grid = linspace(lat_max, lat_min, grid_height); 
lon_grid = linspace(lon_min, lon_max, grid_width);   
[X, Y] = meshgrid(lon_grid, lat_grid); 

%data = readtable('airtelindoordata.csv', 'Delimiter', ',', 'ReadVariableNames', true);


%data = readtable('airtelRoadData.csv', 'Delimiter', ',', 'ReadVariableNames', true);

%data = [data1; data2];
data = readtable('airtelindoordata.csv', 'Delimiter', ',', 'ReadVariableNames', true);
%data = readtable('airtelRoadData.csv', 'Delimiter', ',', 'ReadVariableNames', true);
lon = data.Longitude;
lat = data.Latitude;
rx_values = data.Level; 
valid_idx = ~(isnan(lat) | isnan(lon) | isnan(rx_values));
lon = lon(valid_idx);
lat = lat(valid_idx);
rx_values = rx_values(valid_idx);

x_idx = interp1(lon_grid, 1:grid_width, lon, 'nearest', 'extrap');
y_idx = interp1(lat_grid, 1:grid_height, lat, 'nearest', 'extrap');

valid_points = (x_idx >= 1 & x_idx <= grid_width) & (y_idx >= 1 & y_idx <= grid_height);
lon = lon(valid_points);
lat = lat(valid_points);
rx_values = rx_values(valid_points);


rx_values(rx_values < -140 | rx_values > -30) = NaN;  
valid = ~isnan(rx_values);

F = scatteredInterpolant(lon(valid), lat(valid), rx_values(valid), 'natural', 'none');

rxPowerSmooth = F(X, Y);

rxPowerSmooth = fillmissing(rxPowerSmooth, 'nearest');


figure;

surf(X, Y, rxPowerSmooth, 'EdgeColor', 'none');
colormap jet;
colorbar;
title('Smooth Rx Power Heatmap');
xlabel('Longitude');
ylabel('Latitude');
zlabel('Rx Power (dBm)');
view(2); 

threshold = prctile(rxPowerSmooth(:), 90);

high_power_mask = rxPowerSmooth >= threshold;

CC = bwconncomp(high_power_mask);

stats = regionprops(CC, 'Area', 'Centroid', 'PixelIdxList');

[~, largest_idx] = max([stats.Area]);

centroid_pix = stats(largest_idx).Centroid;

centroid_lon = interp1(1:grid_width, lon_grid, centroid_pix(1));
centroid_lat = interp1(1:grid_height, lat_grid, centroid_pix(2));


region_pixels = stats(largest_idx).PixelIdxList;


region_values = rxPowerSmooth(region_pixels);


[~, max_idx_within_region] = max(region_values);


peak_idx = region_pixels(max_idx_within_region);

[peak_row, peak_col] = ind2sub(size(rxPowerSmooth), peak_idx);

peak_lon = lon_grid(peak_col);
peak_lat = lat_grid(peak_row);

fprintf('Estimated Tower Centroid:\nLongitude: %.6f\nLatitude: %.6f\n', centroid_lon, centroid_lat);
fprintf('Nearest Peak in Red Zone:\nLongitude: %.6f\nLatitude: %.6f\n', peak_lon, peak_lat);


figure;
surf(X, Y, rxPowerSmooth, 'EdgeColor', 'none');
colormap jet;
colorbar;
title('Smooth Rx Power Heatmap with Estimated Tower Location');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Rx Power (dBm)');
view(2);
hold on;


%plot3(centroid_lon, centroid_lat, max(rxPowerSmooth(:)) + 2, 'kp', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
%text(centroid_lon, centroid_lat, max(rxPowerSmooth(:)) + 3, ' Centroid', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

%plot3(peak_lon, peak_lat, max(rxPowerSmooth(:)) + 2, 'mp', 'MarkerSize', 15, 'MarkerFaceColor', 'm');
%text(peak_lon, peak_lat, max(rxPowerSmooth(:)) + 3, ' Peak Tower', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

tower_coordinates.centroid = [centroid_lon, centroid_lat];
tower_coordinates.peak     = [peak_lon, peak_lat];
save('estimated_tower_location_blob_and_peak.mat', 'tower_coordinates');
save('rxPowerSmooth.mat', 'rxPowerSmooth')