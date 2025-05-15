clc 
clear all
grid_width = 1400;
grid_height = 540;

lat_min = 29.861973;  lat_max = 29.864156;
lon_min = 77.895126;  lon_max = 77.901626;

lat_grid = linspace(lat_max, lat_min, grid_height); 
lon_grid = linspace(lon_min, lon_max, grid_width);   
[X, Y] = meshgrid(lon_grid, lat_grid); 


rxPower = -130 * ones(grid_height, grid_width);


% Read the first table
data1 = readtable('airtelindoordata.csv', 'Delimiter', ',', 'ReadVariableNames', true);

% Read another table (e.g. jiooutdoor.csv)
data2 = readtable('airtelRoadData.csv', 'Delimiter', ',', 'ReadVariableNames', true);

data = [data1; data2];
%data = readtable('airtelindoordata.csv', 'Delimiter', ',', 'ReadVariableNames', true);
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
x_idx = x_idx(valid_points);
y_idx = y_idx(valid_points);
rx_values = rx_values(valid_points);

for i = 1:length(x_idx)
    rxPower(y_idx(i), x_idx(i)) = rx_values(i);
end


[Xq, Yq] = meshgrid(1:grid_width, 1:grid_height);
rxPowerSmooth = griddata(x_idx, y_idx, rx_values, Xq, Yq, 'natural');

    
rxPowerSmooth(isnan(rxPowerSmooth)) = -200;


figure(1);

surf(X, Y, rxPower, 'EdgeColor', 'none');
colormap jet;
colorbar;
title('Accurate Rx Power Heatmap');
xlabel('Longitude');
ylabel('Latitude');
zlabel('Rx Power (dBm)');
hold on


figure(2)

opts = detectImportOptions('cellmapperairtel.csv', 'Delimiter', ',', 'ReadVariableNames', false);
cellmap = readtable('cellmapperairtel.csv', opts);

if width(cellmap) == 1
    temp = split(cellmap.Var1, ',');
    cell_lat = str2double(temp(:,1));
    cell_lon = str2double(temp(:,2));
    arfcn    = str2double(temp(:,11));
    tech     = string(temp(:,9));
else
    cell_lat = cellmap.Var1;
    cell_lon = cellmap.Var2;
    arfcn    = cellmap.Var11;
    tech     = string(cellmap.Var9);
end

N = length(cell_lat);


lteBandTable = [
    1,    0,     599,   2110,  1920,   0,      1;
    3,    1200,  1949,  1805,  1710, 1200,    1;
    5,    2400,  2649,  869,   824,  2400,    1;
    40,   38650, 39649, 2300,  2300, 38650,   0;
    41,   39650, 41589, 2496,  NaN,  39650,   0;
];


nrBandTable = [
    78, 620000, 653333, 3300,  0.015, 0;
];


FDL = NaN(N,1);
FUL = NaN(N,1);
Duplex = strings(N,1);
BandNum = NaN(N,1);

for i = 1:N
    if tech(i) == "LTE"
        for b = 1:size(lteBandTable,1)
            N_low  = lteBandTable(b,2);
            N_high = lteBandTable(b,3);
            if arfcn(i) >= N_low && arfcn(i) <= N_high
                FDL(i) = lteBandTable(b,4) + 0.1*(arfcn(i)-lteBandTable(b,6));
                if ~isnan(lteBandTable(b,5))
                    FUL(i) = lteBandTable(b,5) + 0.1*(arfcn(i)-lteBandTable(b,6));
                end
                  if lteBandTable(b,7) == 1
                      Duplex(i) = "FDD";
                 else
                       Duplex(i) = "TDD"; 
                 end
                BandNum(i)= lteBandTable(b,1);
                break;
            end
        end
        
    elseif tech(i) == "NR"
        for b = 1:size(nrBandTable,1)
            N_low  = nrBandTable(b,2);
            N_high = nrBandTable(b,3);
            if arfcn(i) >= N_low && arfcn(i) <= N_high
                FDL(i)    = nrBandTable(b,4) + nrBandTable(b,5)*(arfcn(i)-N_low);
                Duplex(i) = "TDD";
                BandNum(i)= nrBandTable(b,1);
                break;
            end
        end
    end
end


cellInfoTable = table(cell_lat, cell_lon, tech, arfcn, FDL, FUL, Duplex, BandNum);
disp(cellInfoTable);


figure(2);
hold on;

uniqueFDLs = unique(FDL(~isnan(FDL)));
numColors = length(uniqueFDLs);
colors = lines(numColors);

colorIdx = NaN(N,1);
for i = 1:N
    if ~isnan(FDL(i))
        colorIdx(i) = find(FDL(i) == uniqueFDLs);
    end
end

scatterHandles = gobjects(numColors,1);
for c = 1:numColors
    idx = colorIdx == c;
    scatterHandles(c) = scatter(cell_lon(idx), cell_lat(idx), 60, 'x', ...
        'MarkerEdgeColor', colors(c,:), 'LineWidth', 2);
end

legendEntries = strings(numColors,1);
for c = 1:numColors
    legendEntries(c) = sprintf('%.1f MHz (FDL)', uniqueFDLs(c));
end
legend(scatterHandles, legendEntries, 'Location', 'bestoutside');

xlim([lon_min lon_max]);
ylim([lat_min lat_max]);
plot([lon_min lon_max lon_max lon_min lon_min], [lat_min lat_min lat_max lat_max lat_min], 'k--', 'LineWidth', 2);
grid on;
xlabel('Longitude');
ylabel('Latitude');
title('Figure 2: LTE/NR5G ARFCN Mapped Frequency Plot');
hold off;


figure(3);
hold on;

mapImage = imread('map2.png');
imagesc([lon_min lon_max], [lat_min lat_max], flipud(mapImage));
set(gca,'YDir','normal');

scatterHandles = gobjects(numColors,1);
for c = 1:numColors
    idx = colorIdx == c;
    idx_FDD = idx & (Duplex == "FDD");
    idx_TDD = idx & (Duplex == "TDD");

    if any(idx_FDD)
        scatter(cell_lon(idx_FDD), cell_lat(idx_FDD), 60, 'o', ...
            'MarkerEdgeColor', colors(c,:), 'LineWidth', 2);
    end
    if any(idx_TDD)
        scatter(cell_lon(idx_TDD), cell_lat(idx_TDD), 60, '+', ...
            'MarkerEdgeColor', colors(c,:), 'LineWidth', 2);
    end

    scatterHandles(c) = scatter(NaN, NaN, 60, 'x', 'MarkerEdgeColor', colors(c,:), 'LineWidth', 2);
end

legendEntries = strings(numColors,1);
for c = 1:numColors
    legendEntries(c) = sprintf('%.1f MHz (FDL)', uniqueFDLs(c));
end
legend(scatterHandles, legendEntries, 'Location', 'bestoutside');

xlim([lon_min lon_max]);
ylim([lat_min lat_max]);
plot([lon_min lon_max lon_max lon_min lon_min], [lat_min lat_min lat_max lat_max lat_min], 'k--', 'LineWidth', 2);
grid on;
xlabel('Longitude');
ylabel('Latitude');
title('Figure 3: FDL and Duplex Type (o=FDD, +=TDD)');
hold off;


figure(4);
grid on
hold on;

% LTE Bands
numLTEBands = size(lteBandTable,1);
yticksLte = 1:numLTEBands;
yticklabelsLte = strings(numLTEBands,1);

for b = 1:numLTEBands
    
    x_DL = [lteBandTable(b,4), lteBandTable(b,4) + 0.1*(lteBandTable(b,3)-lteBandTable(b,6))];
    fill([x_DL(1), x_DL(2), x_DL(2), x_DL(1)], ...
         [b-0.3, b-0.3, b+0.3, b+0.3], 'yellow', 'EdgeColor', 'k');

    
    if ~isnan(lteBandTable(b,5))
        x_UL = [lteBandTable(b,5), lteBandTable(b,5) + 0.1*(lteBandTable(b,3)-lteBandTable(b,6))];
        fill([x_UL(1), x_UL(2), x_UL(2), x_UL(1)], ...
             [b-0.3, b-0.3, b+0.3, b+0.3], 'cyan', 'EdgeColor', 'k');
    end

   
    duplexType = "TDD";
    if lteBandTable(b,7) == 1
        duplexType = "FDD";
    end
    yticklabelsLte(b) = sprintf('LTE Band %d (%s)', lteBandTable(b,1), duplexType);
end

% NR Bands
numNRBands = size(nrBandTable,1);
yticksNR = numLTEBands+(1:numNRBands);
yticklabelsNR = strings(numNRBands,1);

for b = 1:numNRBands
  
    x_DL = [nrBandTable(b,4), nrBandTable(b,4) + nrBandTable(b,5)*(nrBandTable(b,3)-nrBandTable(b,2))];
    fill([x_DL(1), x_DL(2), x_DL(2), x_DL(1)], ...
         [yticksNR(b)-0.3, yticksNR(b)-0.3, yticksNR(b)+0.3, yticksNR(b)+0.3], ...
         'yellow', 'EdgeColor', 'k');

    
    yticklabelsNR(b) = sprintf('NR Band n%d (TDD)', nrBandTable(b,1));
end


yticks([yticksLte, yticksNR]);
yticklabels([yticklabelsLte; yticklabelsNR]);
xlabel('Frequency (MHz)');
title('Figure 4: Frequency Bands — Downlink (Yellow), Uplink (Cyan)');
grid on;
hold off;

% Step 1: Extract ARFCNs from CellMap table
arfcn_all = cellmap{:, end-1};  % (second last column)

% Step 2: Initialize frequency array
frequencies_all = zeros(size(arfcn_all));

% Step 3: Map each ARFCN to its frequency (FDL)
for idx = 1:length(arfcn_all)
    arfcn = arfcn_all(idx);
    found = false;

    % First check LTE bands
    for b = 1:size(lteBandTable, 1)
        N_low = lteBandTable(b, 2);
        N_high = lteBandTable(b, 3);
        F_low = lteBandTable(b, 4);
        N_off = lteBandTable(b, 6);
        
        if arfcn >= N_low && arfcn <= N_high
            frequencies_all(idx) = F_low + 0.1 * (arfcn - N_off);
            found = true;
            break;
        end
    end

    % If not found, check NR bands
    if ~found
        for b = 1:size(nrBandTable, 1)
            N_low = nrBandTable(b, 2);
            N_high = nrBandTable(b, 3);
            F_low = nrBandTable(b, 4);
            spacing = nrBandTable(b, 5); % in MHz
            % No offset field needed here
            
            if arfcn >= N_low && arfcn <= N_high
                frequencies_all(idx) = F_low + spacing * (arfcn - N_low);
                found = true;
                break;
            end
        end
    end

    % If still not found, set NaN
    if ~found
        frequencies_all(idx) = NaN;
    end
end

% Step 4: Count usage for each uniqueFDL
frequency_counts = zeros(size(uniqueFDLs));
arfcn_labels = strings(size(uniqueFDLs));  % store label text

for i = 1:length(uniqueFDLs)
    % Find matching ARFCNs for this frequency
    idx_matches = find(abs(frequencies_all - uniqueFDLs(i)) < 0.05);
    frequency_counts(i) = numel(idx_matches);
    
    % Get the first matching ARFCN (for labeling)
    if ~isempty(idx_matches)
        arfcn_value = arfcn_all(idx_matches(1));
        arfcn_labels(i) = sprintf('%.1f MHz (ARFCN %d)', uniqueFDLs(i), arfcn_value);
    else
        arfcn_labels(i) = sprintf('%.1f MHz', uniqueFDLs(i));  % fallback
    end
end

% Step 5: Normalize for percentages
frequency_percentages = 100 * frequency_counts / sum(frequency_counts);

% Step 6: Create 3D Pie Chart (Figure 5)
figure(5);
pie3(frequency_percentages);
title('Frequency Usage % in Airtel 5G Data');

% Step 7: Add systematic legend
legend(arfcn_labels, 'Location', 'bestoutside');
