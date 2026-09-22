% rpm_sweep_spin_coat_1d.m
%
% Parameter sweep of spin_coat_stillwagon_topography_solvent_evaporation.mph
% Rotation speed: 400 to 3000 rpm in steps of 200  (14 runs).
%
% Outputs:
%   Figure 1 — h(r) at end time, one curve per rpm
%   Figure 2 — min and mean final film thickness vs rpm
%   Results/rpm_sweep_spin_coat_1d.mat — raw data for all successful runs

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try mphstart; catch; end
import com.comsol.model.util.*

model_path  = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Project\solvent_evaporation_sweep.mph';
results_dir = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Results';
if ~isfolder(results_dir), mkdir(results_dir); end

fprintf('Loading model ...\n');
model = mphload(model_path);
ModelUtil.showProgress(false);

%% Sweep definition
rpm_values = 1000 : 200 : 3000;   % 14 values
n_runs     = length(rpm_values);

%% Storage
r_cell   = cell(1, n_runs);
hh_cell  = cell(1, n_runs);
t_end    = zeros(1, n_runs);
success  = false(1, n_runs);

%% Run loop
for k = 1:n_runs
    rpm = rpm_values(k);
    fprintf('[%2d/%2d]  %4d rpm ... ', k, n_runs, rpm);

    model.param.set('rot', sprintf('%d[rpm]', rpm));

    if k > 1 && success(k-1)
        model.param.set('h_0', sprintf('%.6e[m]', hh_cell{k-1}(end)));
    end

    try
        tic;
        model.study('std1').run();
        elapsed = toc;

        solinfo  = mphsolinfo(model, 'soltag', 'sol1');
        sol_end  = length(solinfo.solvals);
        t_end(k) = solinfo.solvals(sol_end);

        d = mpheval(model, 'hh', 'edim', 1, 'solnum', sol_end);
        r_cell{k}  = d.p(1,:);
        hh_cell{k} = d.d1;
        success(k) = true;

        fprintf('OK   t_end = %5.1f s   h = [%.3f, %.3f] µm   (%.0f s wall)\n', ...
            t_end(k), min(d.d1)*1e6, max(d.d1)*1e6, elapsed);

        save(fullfile(results_dir, 'rpm_sweep_spin_coat_1d.mat'), ...
             'rpm_values', 'r_cell', 'hh_cell', 't_end', 'success');
    catch ME
        fprintf('FAILED — %s\n', ME.message);
    end
end

fprintf('\nCompleted %d / %d runs.\n', sum(success), n_runs);

%% Save results
save(fullfile(results_dir, 'rpm_sweep_spin_coat_1d.mat'), ...
     'rpm_values', 'r_cell', 'hh_cell', 't_end', 'success');
fprintf('Results saved to %s\n', fullfile(results_dir, 'rpm_sweep_spin_coat_1d.mat'));

%% ── Figure 1: h(r) profiles ──────────────────────────────────────────────
figure('Name','RPM Sweep — Film Profile','Position',[50 50 1000 560]);
hold on;

cmap = parula(n_runs);
for k = 1:n_runs
    if success(k)
        lbl = sprintf('%d rpm  (t=%.0f s)', rpm_values(k), t_end(k));
        plot(r_cell{k}*1e3, hh_cell{k}*1e6, ...
             'Color', cmap(k,:), 'LineWidth', 1.8, 'DisplayName', lbl);
    end
end

xlabel('r  [mm]',  'FontSize', 12);
ylabel('h  [µm]',  'FontSize', 12);
title('Film thickness h(r) at end of run — RPM sweep 400–3000 rpm, step 200', ...
      'FontSize', 12);
grid on;
colormap(gca, parula(n_runs));
clim([rpm_values(1), rpm_values(end)]);
cb = colorbar;
cb.Label.String = 'Rotation speed  [rpm]';
cb.Ticks        = linspace(rpm_values(1), rpm_values(end), n_runs);
cb.TickLabels   = arrayfun(@num2str, rpm_values, 'UniformOutput', false);
legend('Location','northeast','FontSize',7,'NumColumns',2);
hold off;

%% ── Figure 2: summary — min and mean h vs rpm ────────────────────────────
idx_ok     = find(success);
rpm_ok     = rpm_values(idx_ok);
hh_min_um  = cellfun(@(v) min(v)*1e6,  hh_cell(idx_ok));
hh_mean_um = cellfun(@(v) mean(v)*1e6, hh_cell(idx_ok));

figure('Name','RPM Sweep — Summary','Position',[50 650 700 380]);
hold on;
plot(rpm_ok, hh_min_um,  'b-o', 'LineWidth', 2, ...
     'MarkerFaceColor', 'b', 'DisplayName', 'min  h');
plot(rpm_ok, hh_mean_um, 'r-s', 'LineWidth', 2, ...
     'MarkerFaceColor', 'r', 'DisplayName', 'mean h');
xlabel('Rotation speed  [rpm]', 'FontSize', 12);
ylabel('Final film thickness  [µm]', 'FontSize', 12);
title('Final film thickness vs rotation speed', 'FontSize', 12);
legend('Location','northeast');
grid on;
hold off;
