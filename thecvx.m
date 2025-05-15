load('final_pathloss_model.mat');  % Loads A, B, peak_lon, peak_lat, rxPowerSmooth
d = linspace(1, 1000, 1000);  
pl_exact = A + B * log10(d);  

p = polyfit(d, log10(d), 1);  
log10_approx = polyval(p, d);
pl_approx = A + B * log10_approx;

figure;
plot(d, pl_exact, 'b', 'LineWidth', 2); hold on;
plot(d, pl_approx, 'r--', 'LineWidth', 2);
legend('Exact: A + B log_{10}(d)', 'Approx: A + B Poly(d)');
xlabel('Distance d (m)');
ylabel('Path Loss (dB)');
title('Path Loss vs Distance: Exact vs Polynomial Approximation');
grid on;

nCols = 1400;
nRows = 540;

lat_max = 29.8642077;  
lat_min = 29.861973;   
lon_min = 77.895126;   
lon_max = 77.901626;  

[x_grid, y_grid] = meshgrid(0:nCols-1, 0:nRows-1);

lon_grid = lon_min + (x_grid / (nCols - 1)) * (lon_max - lon_min);
lat_grid = lat_min + (y_grid / (nRows - 1)) * (lat_max - lat_min);  % Flip y

x_flat = x_grid(:);
y_flat = y_grid(:);
n = length(x_flat);

x1 = round((peak_lon - lon_min) / (lon_max - lon_min) * (nCols - 1));
y1 = round((peak_lat - lat_min) / (lat_max - lat_min) * (nRows - 1));

center_lat = (lat_min + lat_max) / 2;
dlat_m = deg2rad(lat_max - lat_min) * 6371000;
dlon_m = deg2rad(lon_max - lon_min) * 6371000 * cosd(center_lat);
dist_per_pix_y = dlat_m / (nRows - 1);
dist_per_pix_x = dlon_m / (nCols - 1);

figure;
plot(x_flat, y_flat, 'k.', 'MarkerSize', 1); hold on;
plot(x1, y1, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(x1+10, y1, 'Tower 1', 'Color', 'r', 'FontSize', 10);
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Meshgrid and Tower 1 Location'); grid on; axis equal;
tx_power = 53;  % Transmit power (adjust as needed)

cvx_begin
cvx_solver mosek
    variables x2 dx1 dy1 dx2 dy2 total_rx_power d path_loss1 path_loss2 d_meters1 d_meters2
    variable y2

    % Distance from Tower 1 to all pixels
    dx1 >= (x_flat - x1) * dist_per_pix_x;
    dy1 >= (y_flat - y1) * dist_per_pix_y;
    d_meters1 >= norm([dx1, dy1]);  % Distance to Tower 1

    % Distance from Tower 2 to all pixels
    dx2 >= (x_flat - x2) * dist_per_pix_x;
    dy2 >= (y_flat - y2) * dist_per_pix_y;
    d_meters2 >= norm([dx2, dy2]);  % Distance to Tower 2
    
    d <= 1500
    
    d >= norm([x2-x1, y2-y1])

    % Polynomial approximation (degree 1) for log10(d)
    log10_approx1 = p(1) * d_meters1 + p(2);  % Applying the polynomial for Tower 1
    log10_approx2 = p(1) * d_meters2 + p(2);  % Applying the polynomial for Tower 2
    600 <= x2 <= nCols-1;
    0 <= y2 <= nRows-1;
    % Path loss models using polynomial approximation
    path_loss1 >= A + B * log10_approx1;  % Path loss for Tower 1
    path_loss2 >= A + B * log10_approx2;  % Path loss for Tower 2

    % Total received power from both towers
    total_rx_power >= max((tx_power - path_loss1) ,  (tx_power - path_loss2));

    % Threshold condition
    rx_power_threshold = -90;

    % Minimize the total received power to meet the threshold
    minimize(total_rx_power - 2*d)
    
    subject to
        total_rx_power >= rx_power_threshold +10 ;  % Total power must be above threshold
        path_loss1 + path_loss2 >= 0
        
        path_loss1 >= 0
        path_loss2 >= 0
        d >= 0
        
cvx_end


lon2 = lon_min + (x2 / (nCols - 1)) * (lon_max - lon_min);
lat2 = lat_max - (y2 / (nRows - 1)) * (lat_max - lat_min);
tx2 = tx_power;

R = 6371000; 
dlat = deg2rad(lat_grid - lat2);
dlon = deg2rad(lon_grid - lon2);
a = sin(dlat/2).^2 + cosd(lat_grid) .* cosd(lat2) .* sin(dlon/2).^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));
d2 = R * c;

PL2 = A + B * log10(d2 + 1);
rx2 = tx2 - PL2;
final_rx = max(rxPowerSmooth, rx2);

figure;
surf(lon_grid, lat_grid, final_rx, 'EdgeColor', 'none');
view(2); colormap jet; colorbar;
title('Final Combined Rx Power with Two Towers');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Rx Power (dBm)'); hold on;
plot3(peak_lon, peak_lat, max(final_rx(:)) + 1, 'mp', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
text(peak_lon, peak_lat, max(final_rx(:)) + 2, ' Tower 1', 'Color', 'w', 'FontSize', 9);
plot3(lon2, lat2, max(final_rx(:)) + 1, 'gp', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
text(lon2, lat2, max(final_rx(:)) + 2, ' Tower 2', 'Color', 'w', 'FontSize', 9);
caxis([-110 -50]);  % <-- Fix color scale from -120 dBm (blue) to -40 dBm (red)

hold off;

save('second_tower_results.mat', 'lon2', 'lat2', 'tx2')