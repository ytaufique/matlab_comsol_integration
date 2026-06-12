% Parameters
x1 = 0.47;
x2 = 0.50;
y1 = 0.50;
y2 = 0.53;
A  = 10e-6;
epsilon = 3e-4;

% Grid
[x, y] = meshgrid(linspace(0, 1, 500), linspace(0, 1, 500));

flc2hs = @(u, e) (u < -e) .* 0 + ...
                 (u >  e) .* 1 + ...
                 (abs(u) <= e) .* (3*(u/e + 1).^2/4 - (u/e + 1).^3/4);

% Square bump
s = A * (flc2hs(x - x1, epsilon) - flc2hs(x - x2, epsilon)) .* ...
        (flc2hs(y - y1, epsilon) - flc2hs(y - y2, epsilon));

% Plot
figure;
surf(x, y, s, 'EdgeColor', 'none');
colorbar;
xlabel('x'); ylabel('y'); zlabel('s');
title('Square Bump');
view(3);