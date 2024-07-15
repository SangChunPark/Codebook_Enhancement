grid_theta  = 90 - 90*(1 + margin) : gap : 90 + 90*(1 + margin);
grid_phi    = -90*(1 + margin) : gap : 90*(1 + margin);

[x1_, x2_] = meshgrid(grid_theta, grid_phi);
x1 = x1_(:);
x2 = x2_(:);
xi = [x1 x2];
xj = xi;

xj(xj(:, 1) < 0, 1)     = xj(xj(:, 1) < 0, 1) + 180;
xj(xj(:, 1) >= 180, 1)  = xj(xj(:, 1) >= 180, 1) - 180;

xj(xj(:, 2) < -90, 2)   = xj(xj(:, 2) < -90, 2) + 180;
xj(xj(:, 2) >= 90, 2)   = xj(xj(:, 2) >= 90, 2) - 180;

xj = round(xj, 3);
[xj, ~, ic] = unique(xj, 'rows');