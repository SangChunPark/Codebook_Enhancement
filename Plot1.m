figure(iteration + 1)
histogram2(x(:, 1), x(:, 2), 'BinWidth', [3 3])
xlim([min(grid_theta) max(grid_theta)])
ylim([min(grid_phi) max(grid_phi)])

xlabel('Zenith angle [deg]')
ylabel('Azimuth angle [deg]')
zlabel('The number of PMI indices')
xticks(0 : 30 : 180)
yticks(-90 : 30 : 90)
xlim([0 180])
ylim([-90 90])
set(gca, 'PlotBoxAspectRatio', [1, 1, 0.5])
figure(iteration + 21)
x = unique(x, 'rows');

title(['Directions of PMI indices in iteration ', num2str(iteration)])
plot(x(:, 2), 180 - x(:, 1), 's', 'markersize', 5, 'linewidth', 1.5)
xlabel('Azimuth angle [deg]')
ylabel('Zenith angle [deg]')
xlim([-90 90])
ylim([0 180])
xticks(-90 : 30 : 90)
yticks(0 : 30 : 180)
yticklabels({'180', '150', '120', '90', '60', '30', '0'})
grid
set(gca, 'PlotBoxAspectRatio', [1, 1, 1])