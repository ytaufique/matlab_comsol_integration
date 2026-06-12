% run_spin_coat_2d.m
%
% 1. Builds the 2D spin-coat model (runs create_spin_coat_2d.m if the .mph
%    does not exist yet, or reloads it if it does).
% 2. Solves the transient study (t = 0..100 s).
% 3. Extracts hh, P, c on a regular grid at t = 99 s.
% 4. Produces a 3D surface visualisation of film thickness.

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try mphstart; catch; end
import com.comsol.model.util.*

model_path = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Project\spin_coat_2d.mph';
script_dir = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Matlab';

%% ── 1. Build or load the model ────────────────────────────────────────────
if isfile(model_path)
    fprintf('Loading existing model: %s\n', model_path);
    model = mphload(model_path);
else
    fprintf('Model not found — running create_spin_coat_2d.m ...\n');
    run(fullfile(script_dir, 'create_spin_coat_2d.m'));
    fprintf('Model built and saved.\n');
end

%% ── 2. Solve ──────────────────────────────────────────────────────────────
ModelUtil.showProgress(true);
fprintf('\nSolving ... (t = 0 to 100 s)\n');
tic;
model.study('std1').run();
fprintf('Solve complete (%.1f s wall time).\n\n', toc);

%% ── 3. Extract results at t = 99 s on a regular 2D grid ──────────────────
L    = str2double(char(model.param.get('L')));   % domain half-size [m]
t_plot = 99;                                      % extraction time [s]

Nplot = 150;                                      % grid points per axis
xv    = linspace(-L, L, Nplot);
yv    = linspace(-L, L, Nplot);
[Xp, Yp] = meshgrid(xv, yv);
coords    = [Xp(:)'; Yp(:)'];                     % 2 x N_points

fprintf('Interpolating hh, P, c at t = %g s on %dx%d grid ...\n', ...
        t_plot, Nplot, Nplot);

% mphinterp interpolates any expression at arbitrary (x,y) coordinates
hh_vec = mphinterp(model, 'hh', 'coord', coords, ...
                   'dataset', 'dset1', 't', t_plot);
P_vec  = mphinterp(model, 'P',  'coord', coords, ...
                   'dataset', 'dset1', 't', t_plot);
c_vec  = mphinterp(model, 'c',  'coord', coords, ...
                   'dataset', 'dset1', 't', t_plot);

HH = reshape(hh_vec, Nplot, Nplot) * 1e6;   % [m] -> [µm]
PP = reshape(P_vec,  Nplot, Nplot);          % [Pa]
CC = reshape(c_vec,  Nplot, Nplot);          % [-]

xmm = xv * 1e3;   % [m] -> [mm]
ymm = yv * 1e3;

fprintf('  hh range: %.3f – %.3f µm\n', min(HH(:)), max(HH(:)));
fprintf('  P  range: %.3f – %.3f Pa\n',  min(PP(:)), max(PP(:)));
fprintf('  c  range: %.4f – %.4f\n',     min(CC(:)), max(CC(:)));

%% ── 4. 3D surface: film thickness h(x,y) at t = 99 s ────────────────────
figure('Name','2D Spin Coat — Film Thickness', ...
       'Position', [60 60 1050 780]);

%  ── Main 3D surface ──────────────────────────────────────────────────────
ax1 = subplot(2,2,[1 3]);
surf(xmm, ymm, HH, 'EdgeColor', 'none');
colormap(ax1, 'parula');
cb = colorbar;
cb.Label.String = 'h  [\mum]';
xlabel('x  [mm]');
ylabel('y  [mm]');
zlabel('h  [\mum]');
title(sprintf('Film thickness  h(x,y,t=%.0f s)', t_plot), 'FontSize', 12);
view(-35, 30);
axis tight;
set(ax1, 'FontSize', 10);

% Overlay the topography ring for reference
hold on;
offset_mm = str2double(char(model.param.get('offset'))) * 1e3;
theta_ring = linspace(0, 2*pi, 200);
z_ring     = ones(1,200) * min(HH(:)) - 0.01*(max(HH(:))-min(HH(:)));
plot3(offset_mm*cos(theta_ring), offset_mm*sin(theta_ring), z_ring, ...
      'w--', 'LineWidth', 1.2, 'DisplayName', 'Topography ring');
hold off;
legend('Location','northeast');

%  ── Top-view (filled contour) ────────────────────────────────────────────
ax2 = subplot(2,2,2);
contourf(xmm, ymm, HH, 30, 'LineColor', 'none');
colormap(ax2, 'parula');
colorbar;
axis equal tight;
xlabel('x  [mm]');  ylabel('y  [mm]');
title('Top view (filled contour)', 'FontSize', 10);
hold on;
plot(offset_mm*cos(theta_ring), offset_mm*sin(theta_ring), ...
     'w--', 'LineWidth', 1.0);
hold off;

%  ── Radial profile along y=0 ─────────────────────────────────────────────
ax3 = subplot(2,2,4);
[~, iy0] = min(abs(yv));
plot(xmm, HH(iy0,:), 'b-', 'LineWidth', 2);
hold on;
xline(offset_mm,  'k--', 'Topog.', 'LabelHorizontalAlignment','left');
xline(-offset_mm, 'k--');
hold off;
xlabel('x  [mm]');
ylabel('h  [\mum]');
title('Radial profile  h(x, y=0)', 'FontSize', 10);
grid on;

sgtitle(sprintf('2D Spin Coat Simulation  |  t = %.0f s', t_plot), ...
        'FontSize', 13, 'FontWeight', 'bold');

%% ── 5. Supplementary: concentration c(x,y) ──────────────────────────────
figure('Name','2D Spin Coat — Concentration & Pressure', ...
       'Position', [1130 60 900 420]);

subplot(1,2,1);
surf(xmm, ymm, CC, 'EdgeColor','none');
colormap('turbo');  colorbar;
xlabel('x [mm]');  ylabel('y [mm]');  zlabel('c [-]');
title(sprintf('Concentration  c(x,y,t=%.0f s)', t_plot));
view(-35, 30);  axis tight;

subplot(1,2,2);
surf(xmm, ymm, PP, 'EdgeColor','none');
colormap('coolwarm');  colorbar;
xlabel('x [mm]');  ylabel('y [mm]');  zlabel('P [Pa]');
title(sprintf('Pressure  P(x,y,t=%.0f s)', t_plot));
view(-35, 30);  axis tight;

sgtitle(sprintf('t = %.0f s', t_plot), 'FontSize', 12);
