clc
clear all
d = linspace(0.0001, 500,500);


PL_indoor = 105.97 + 13.19 * log10(d);

PL_3gppInhLos = 43 + 17.3 * log10(d); 

figure;
plot(d, PL_3gppInhLos, 'r--', 'LineWidth', 2); hold on;
plot(d, PL_indoor, 'b', 'LineWidth', 2);
legend('InH - Office LOS Model', 'Observed Indoor Model', 'Location', 'northwest');
xlabel('Distance d (m)');
ylabel('Path Loss (dB)');
title('Indoor: Observed vs 3GPP Path Loss');
grid on;

PL_road = 79.79 + 22.85 * log10(d);
PL_3gppUmaLOS = 38.6 + 22 * log10(d);  

figure;
plot(d, PL_3gppUmaLOS, 'r--', 'LineWidth', 2); hold on;
plot(d, PL_road, 'b', 'LineWidth', 2);
legend('UMa LOS', 'Observed Road Model', 'Location', 'northwest');
xlabel('Distance d (m)');
ylabel('Path Loss (dB)');
title('Road: Observed vs 3GPP Path Loss');
grid on;
