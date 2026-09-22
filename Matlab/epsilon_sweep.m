% epsilon_sweep.m
%
% Parameter sweep of solvent_evaporation_sweep_solvent.mph over epsilon [m].
% Runs continue even if individual simulations crash.
%
% Outputs:
%   Results/epsilon_sweep.mat — r_cell, hh_cell, s_cell, t_end, success
%   Figure 1 — max slope of (hh+s) vs max slope of s, one point per epsilon

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try mphstart; catch; end
import com.comsol.model.util.*

model_path  = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Project\solvent_evaporation_sweep_solvent.mph';
results_dir = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Results';
if ~isfolder(results_dir), mkdir(results_dir); end

fprintf('Loading model ...\n');
model = mphload(model_path);
ModelUtil.showProgress(false);

%% Sweep definition
% s = A * flc2hs(x_prime, epsilon), A = 2e-6 m
% max slope of s = atan(A*15/(16*epsilon)) * 180/pi
% 5 deg -> epsilon ~ 21.4e-6 m,  20 deg -> epsilon ~ 5.1e-6 m
epsilon_values = linspace(5.1e-6, 21.4e-6, 7);   % [m]
n_runs = length(epsilon_values);

%% Storage
r_cell  = cell(1, n_runs);
hh_cell = cell(1, n_runs);
s_cell  = cell(1, n_runs);
t_end   = zeros(1, n_runs);
success = false(1, n_runs);

%% Run loop
for k = 1:n_runs
    epsilon = epsilon_values(k);
    fprintf('[%2d/%2d]  epsilon = %.2e m ... ', k, n_runs, epsilon);

    model.param.set('epsilon', sprintf('%.6e[m]', epsilon));

    try
        tic;
        model.study('std1').run();
        elapsed = toc;

        solinfo    = mphsolinfo(model, 'soltag', 'sol1');
        sol_end    = length(solinfo.solvals);
        t_end(k)   = solinfo.solvals(sol_end);

        d           = mpheval(model, 'hh', 'edim', 1, 'solnum', sol_end);
        r_cell{k}  = d.p(1,:);
        hh_cell{k} = d.d1;

        ds          = mpheval(model, 's',  'edim', 1, 'solnum', sol_end);
        s_cell{k}  = ds.d1;

        success(k) = true;

        fprintf('OK   t_end = %5.1f s   h = [%.3f, %.3f] µm   (%.0f s wall)\n', ...
            t_end(k), min(d.d1)*1e6, max(d.d1)*1e6, elapsed);

        save(fullfile(results_dir, 'epsilon_sweep.mat'), ...
             'epsilon_values', 'r_cell', 'hh_cell', 's_cell', 't_end', 'success');
    catch ME
        fprintf('FAILED — %s\n', ME.message);
    end
end

fprintf('\nCompleted %d / %d runs.\n', sum(success), n_runs);

%% ── Figure 1: max slope of (hh+s) vs max slope of s ─────────────────────
idx_ok        = find(success);
max_slope_s   = zeros(1, numel(idx_ok));
max_slope_hhs = zeros(1, numel(idx_ok));

for j = 1:numel(idx_ok)
    k  = idx_ok(j);
    dr = diff(r_cell{k});

    slope_s        = atan(diff(s_cell{k})  ./ dr) * 180 / pi;
    max_slope_s(j) = max(abs(slope_s));

    slope_hhs         = atan(diff(hh_cell{k} + s_cell{k}) ./ dr) * 180 / pi;
    max_slope_hhs(j)  = max(abs(slope_hhs));
end

figure('Name','Epsilon Sweep — Max slope','Position',[50 50 700 480]);
scatter(max_slope_s, max_slope_hhs, 80, log10(epsilon_values(idx_ok)), 'filled');
cb = colorbar;
cb.Label.String = 'log_{10}(epsilon  [m])';
xlabel('Max slope of s  [°]',     'FontSize', 12);
ylabel('Max slope of hh+s  [°]',  'FontSize', 12);
title('Film surface slope vs substrate slope — \epsilon sweep', 'FontSize', 12);
grid on;

for j = 1:numel(idx_ok)
    text(max_slope_s(j), max_slope_hhs(j), ...
         sprintf('  %.0e m', epsilon_values(idx_ok(j))), ...
         'FontSize', 8, 'VerticalAlignment', 'middle');
end
saveas(gcf, fullfile(results_dir, 'epsilon_sweep_max_slope.png'));
