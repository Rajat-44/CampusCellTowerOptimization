% === Load tower 2 result and model coefficients ===
load('final_pathloss_model.mat');  % A, B, peak_lon, peak_lat, rxPowerSmooth
load('second_tower_results.mat');  % lat2, lon2, tx2

nCols = 1400;
nRows = 540;

lat_min = 29.861973;  lat_max = 29.8642077;
lon_min = 77.895126;  lon_max = 77.901626;

[x_grid, y_grid] = meshgrid(0:nCols-1, 0:nRows-1);
lon_grid = lon_min + (x_grid / (nCols - 1)) * (lon_max - lon_min);
lat_grid = lat_max - (y_grid / (nRows - 1)) * (lat_max - lat_min);  % flip vertically

R = 6371000;
dlat = deg2rad(lat_grid - lat2);
dlon = deg2rad(lon_grid - lon2);
a = sin(dlat/2).^2 + cosd(lat_grid) .* cosd(lat2) .* sin(dlon/2).^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));
d2 = R * c;  % meters

PL2 = A + B * log10(d2);
rx2 = tx2 - PL2;

final_rx = max(rxPowerSmooth, rx2);

figure;
surf(lon_grid, lat_grid, rxPowerSmooth, 'EdgeColor', 'none');
colormap jet; colorbar;
title('Rx Power from Tower 1');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Rx Power (dBm)');
caxis([-110 -50]);  % <-- Fix color scale from -120 dBm (blue) to -40 dBm (red)
view(2);

figure;
surf(lon_grid, lat_grid, rx2, 'EdgeColor', 'none');
colormap jet; colorbar;
title('Rx Power from Tower 2 (Direct Log Model)');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Rx Power (dBm)');
caxis([-110 -50]);  % <-- Fix color scale from -120 dBm (blue) to -40 dBm (red)

view(2);

figure;
surf(lon_grid, lat_grid, final_rx, 'EdgeColor', 'none');
colormap jet; colorbar;
title('Maximum Rx Power with Two Towers');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Rx Power (dBm)');
caxis([-110 -50]);  % <-- Fix color scale from -120 dBm (blue) to -40 dBm (red)

view(2); hold on;

plot3(peak_lon, peak_lat, max(final_rx(:)) + 1, 'mp', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
text(peak_lon, peak_lat, max(final_rx(:)) + 2, ' Tower 1', 'Color', 'w', 'FontSize', 9);
plot3(lon2, lat2, max(final_rx(:)) + 1, 'gp', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
text(lon2, lat2, max(final_rx(:)) + 2, ' Tower 2', 'Color', 'w', 'FontSize', 9);
hold off;

save('rxPowerTower2_logbased.mat', 'rx2');
save('combined_rxPower_logbased.mat', 'final_rx', 'lon_grid', 'lat_grid');


hold on;
plot3(peak_lon, peak_lat, max(final_rx(:)) + 1, 'mp', 'MarkerSize', 12, 'MarkerFaceColor', 'm');  % Tower 1
text(peak_lon, peak_lat, max(final_rx(:)) + 2, ' Tower 1', 'Color', 'w', 'FontSize', 9);

plot3(lon2, lat2, max(final_rx(:)) + 1, 'gp', 'MarkerSize', 12, 'MarkerFaceColor', 'g');  % Tower 2
text(lon2, lat2, max(final_rx(:)) + 2, ' Tower 2', 'Color', 'w', 'FontSize', 9);
hold off;
mapImage = imread('map2.png'); 
[imgH, imgW, ~] = size(mapImage);
fprintf('Image size: %d x %d\n', imgW, imgH);

figure;
imagesc([lon_min lon_max], [lat_min lat_max], flipud(mapImage)); 
set(gca,'YDir','normal');
xlabel('Longitude'); ylabel('Latitude');
title('Tower Locations on Map');
hold on;

plot(peak_lon, peak_lat, 'mp', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
text(peak_lon + 0.00015, peak_lat, 'Tower 1', ...
     'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold', ...
     'BackgroundColor', 'w', 'Margin', 2, 'EdgeColor', 'k', 'Clipping', 'on');

plot(lon2, lat2, 'gp', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
text(lon2 + 0.00015, lat2, 'Tower 2', ...
     'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold', ...
     'BackgroundColor', 'w', 'Margin', 2, 'EdgeColor', 'k', 'Clipping', 'on');

hold off;