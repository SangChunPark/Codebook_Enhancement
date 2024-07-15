figure(iteration + 61)

plot(enhanced.sph(2, :), 180 - enhanced.sph(1, :), 's', 'markersize', 5, 'linewidth', 1.5)
xlabel('Azimuth angle [deg]')
ylabel('Zenith angle [deg]')
xlim([-90 90])
ylim([0 180])
xticks(-90 : 30 : 90)
yticks(0 : 30 : 180)
yticklabels({'180', '150', '120', '90', '60', '30', '0'})
grid
set(gca, 'PlotBoxAspectRatio', [1, 1, 1])