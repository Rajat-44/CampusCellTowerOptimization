clc; clear; close all;

 
grid_width = 1400;
grid_height = 540;

lat_min = 29.861973;  lat_max = 29.864156;
lon_min = 77.895126;  lon_max = 77.901626;

lat_grid = linspace(lat_max, lat_min, grid_height);
lon_grid = linspace(lon_min, lon_max, grid_width);
[X, Y] = meshgrid(lon_grid, lat_grid);


%data = readtable('jioindoordata.csv', 'Delimiter', ',', 'ReadVariableNames', true);
data = readtable('airtelindoordata.csv', 'Delimiter', ',', 'ReadVariableNames', true);
%data = readtable('jioRoadData.csv', 'Delimiter', ',', 'ReadVariableNames', true);
%data = readtable('airtelRoadData.csv', 'Delimiter', ',', 'ReadVariableNames', true);

%data = [data1; data2];

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


threshold = prctile(rxPowerSmooth(:), 90);
high_power_mask = rxPowerSmooth >= threshold;
CC = bwconncomp(high_power_mask);
stats = regionprops(CC, 'Area', 'Centroid', 'PixelIdxList');
[~, largest_idx] = max([stats.Area]);
region_pixels = stats(largest_idx).PixelIdxList;
region_values = rxPowerSmooth(region_pixels);
[~, max_idx_within_region] = max(region_values);
peak_idx = region_pixels(max_idx_within_region);
[peak_row, peak_col] = ind2sub(size(rxPowerSmooth), peak_idx);
peak_lon = lon_grid(peak_col);
peak_lat = lat_grid(peak_row);

% jio tower 29.862661 77.896217
peak_lon = 77.896217;
peak_lat = 29.862661;
%
%airtel
%

fprintf('Estimated Tower Location:\nLongitude: %.6f\nLatitude: %.6f\n', peak_lon, peak_lat);

R = 6371000; 
dlat = deg2rad(lat_grid' - peak_lat);
dlon = deg2rad(lon_grid - peak_lon);
a = sin(dlat/2).^2 + cos(deg2rad(peak_lat)) * cos(deg2rad(lat_grid')) .* sin(dlon/2).^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
D = R * c;

Pt = 43; % Assumed Tx power dBm
PL_grid = Pt - rxPowerSmooth;

dist_vec = D(:);
pl_vec = PL_grid(:);
valid_idx = isfinite(pl_vec) & isfinite(dist_vec) & dist_vec > 1;
dist_vec_clean = dist_vec(valid_idx);
pl_vec_clean = pl_vec(valid_idx);

modelfun = @(b, d) b(1) + b(2)*log10(d);
b0 = [30, 20];
mdl = fitnlm(dist_vec_clean, pl_vec_clean, modelfun, b0);

A = mdl.Coefficients.Estimate(1);
B = mdl.Coefficients.Estimate(2);

fprintf('\nFitted Path Loss Model:\n');
fprintf('PL(d) = %.2f + %.2f * log10(d)\n', A, B);

dist_range = linspace(min(dist_vec_clean), max(dist_vec_clean), 200);
pl_pred = A + B*log10(dist_range);

figure;
scatter(dist_vec_clean, pl_vec_clean, 8, 'b', 'filled'); hold on;
plot(dist_range, pl_pred, 'r-', 'LineWidth', 2);
grid on;
xlabel('Distance from Tower (m)');
ylabel('Path Loss (dB)');
title('Measured Path Loss and Fitted Log-Distance Model');
legend('Measured Path Loss', 'Fitted Model');
xlim([0 max(dist_vec_clean)]);
ylim([min(pl_vec_clean)-5, max(pl_vec_clean)+5]);

%% Plot Heatmap with Tower Location and 261m, 300m Circles
figure;
surf(X, Y, rxPowerSmooth, 'EdgeColor', 'none');
colormap jet; colorbar;
title('Smoothed Rx Power Heatmap with Tower Location + 261m & 300m Circles');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Rx Power (dBm)');
view(2);
hold on;

% Mark tower peak
plot3(peak_lon, peak_lat, max(rxPowerSmooth(:))+2, 'mp', 'MarkerSize', 15, 'MarkerFaceColor', 'm');
text(peak_lon, peak_lat, max(rxPowerSmooth(:))+3, ' Peak Tower', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

% Earth radius and theta for circles
R_earth = 6371000; % meters
theta = linspace(0, 2*pi, 300);

% ===== 261 m Circle =====
d_261 = 261;
delta_lat_261 = (d_261 / R_earth) * (180/pi);
delta_lon_261 = delta_lat_261 / cosd(peak_lat);
circle_lat_261 = peak_lat + (delta_lat_261 * cos(theta));
circle_lon_261 = peak_lon + (delta_lon_261 * sin(theta));

plot3(circle_lon_261, circle_lat_261, repmat(max(rxPowerSmooth(:))+1, size(circle_lon_261)), 'w-', 'LineWidth', 2);
text(peak_lon, peak_lat + delta_lat_261, max(rxPowerSmooth(:))+1.5, ' 261m', 'Color', 'b', 'FontSize', 9, 'FontWeight', 'bold');

% ===== 300 m Circle =====
d_300 = 300;
delta_lat_300 = (d_300 / R_earth) * (180/pi);
delta_lon_300 = delta_lat_300 / cosd(peak_lat);
circle_lat_300 = peak_lat + (delta_lat_300 * cos(theta));
circle_lon_300 = peak_lon + (delta_lon_300 * sin(theta));

plot3(circle_lon_300, circle_lat_300, repmat(max(rxPowerSmooth(:))+1, size(circle_lon_300)), 'w-', 'LineWidth', 2);
text(peak_lon, peak_lat + delta_lat_300, max(rxPowerSmooth(:))+1.5, ' 300m', 'Color', 'b', 'FontSize', 9, 'FontWeight', 'bold');

hold off;

save('final_pathloss_model.mat', 'mdl', 'A', 'B', 'peak_lon', 'peak_lat', 'rxPowerSmooth');

