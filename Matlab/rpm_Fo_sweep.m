% rpm_Fo_sweep.m
%
% 2-D parameter sweep of solvent_evaporation_sweep_solvent.mph
% Outer loop: evaporation flux Fo  [m/s]
% Inner loop: rotation speed rot   [rpm]
%
% Outputs:
%   Results/rpm_Fo_sweep.mat — r_cell, hh_cell, t_end, success (all n_Fo x n_rpm)

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
Fo_values  = [2e-8, 5e-8, 10e-8, 20e-8];   % evaporation flux [m/s]
rpm_values = 2000 : 200 : 3000;                     % rotation speed [rpm]
n_Fo  = length(Fo_values);
n_rpm = length(rpm_values);
n_total = n_Fo * n_rpm;

%% Storage  (rows = Fo, columns = rpm)
r_cell  = cell(n_Fo, n_rpm);
hh_cell = cell(n_Fo, n_rpm);
s_cell  = cell(n_Fo, n_rpm);
t_end   = zeros(n_Fo, n_rpm);
success = false(n_Fo, n_rpm);

%% Run loop
run_idx = 0;
for i = 1:n_Fo
    Fo  = Fo_values(i);
    model.param.set('Fo', sprintf('%.6e[m/s]', Fo));

    for k = 1:n_rpm
        rpm = rpm_values(k);
        run_idx = run_idx + 1;
        fprintf('[%3d/%3d]  Fo = %.1e m/s   %4d rpm ... ', run_idx, n_total, Fo, rpm);

        model.param.set('rot', sprintf('%d[rpm]', rpm));

        try
            tic;
            model.study('std1').run();
            elapsed = toc;

            solinfo       = mphsolinfo(model, 'soltag', 'sol1');
            sol_end       = length(solinfo.solvals);
            t_end(i,k)    = solinfo.solvals(sol_end);

            d              = mpheval(model, 'hh', 'edim', 1, 'solnum', sol_end);
            r_cell{i,k}   = d.p(1,:);
            hh_cell{i,k}  = d.d1;

            ds             = mpheval(model, 's',  'edim', 1, 'solnum', sol_end);
            s_cell{i,k}   = ds.d1;

            success(i,k)  = true;

            fprintf('OK   t_end = %5.1f s   h = [%.3f, %.3f] µm   (%.0f s wall)\n', ...
                t_end(i,k), min(d.d1)*1e6, max(d.d1)*1e6, elapsed);

            save(fullfile(results_dir, 'rpm_Fo_sweep.mat'), ...
                 'Fo_values', 'rpm_values', 'r_cell', 'hh_cell', 's_cell', 't_end', 'success');
        catch ME
            fprintf('FAILED — %s\n', ME.message);
        end
    end
end

fprintf('\nCompleted %d / %d runs.\n', sum(success(:)), n_total);

%% ── Figure 1: mean h vs rpm, one curve per Fo ───────────────────────────
figure('Name','RPM-Fo Sweep — Mean h vs rpm','Position',[50 50 800 500]);
hold on;
cmap = parula(n_Fo);
for i = 1:n_Fo
    idx_ok = find(success(i,:));
    if isempty(idx_ok), continue; end
    hh_mean = cellfun(@(v) mean(v)*1e6, hh_cell(i, idx_ok));
    plot(rpm_values(idx_ok), hh_mean, '-o', 'Color', cmap(i,:), ...
         'LineWidth', 2, 'MarkerFaceColor', cmap(i,:), ...
         'DisplayName', sprintf('Fo = %.0e m/s', Fo_values(i)));
end
xlabel('Rotation speed  [rpm]', 'FontSize', 12);
ylabel('Mean final film thickness  [µm]', 'FontSize', 12);
title('Mean h vs rotation speed — Fo sweep', 'FontSize', 12);
legend('Location', 'northeast');
grid on;
hold off;

%% ── Figure 2: heatmap — mean h over Fo x rpm grid ───────────────────────
hh_mean_grid = zeros(n_Fo, n_rpm);
for i = 1:n_Fo
    for k = 1:n_rpm
        if success(i,k)
            hh_mean_grid(i,k) = mean(hh_cell{i,k}) * 1e6;
        else
            hh_mean_grid(i,k) = NaN;
        end
    end
end

figure('Name','RPM-Fo Sweep — Heatmap','Position',[900 50 700 450]);
imagesc(rpm_values, 1:n_Fo, hh_mean_grid);
set(gca, 'YTick', 1:n_Fo, 'YTickLabel', arrayfun(@(v) sprintf('%.0e', v), Fo_values, 'UniformOutput', false));
xlabel('Rotation speed  [rpm]', 'FontSize', 12);
ylabel('Fo  [m/s]', 'FontSize', 12);
title('Mean final film thickness  [µm]', 'FontSize', 12);
colorbar;
colormap(parula);

%% ── Figure 3: max hh slope vs Fo, one line per rpm ──────────────────────
max_slope = NaN(n_Fo, n_rpm);
for i = 1:n_Fo
    for k = 1:n_rpm
        if success(i,k)
            dr    = diff(r_cell{i,k});
            dhh   = diff(hh_cell{i,k} + s_cell{i,k});
            slope = atan(dhh ./ dr) * 180 / pi;
            max_slope(i,k) = max(abs(slope));
        end
    end
end

figure('Name','RPM-Fo Sweep — Max slope vs Fo','Position',[50 550 800 450]);
hold on;
cmap2 = parula(n_rpm);
for k = 1:n_rpm
    idx_ok = find(~isnan(max_slope(:,k)));
    if isempty(idx_ok), continue; end
    plot(Fo_values(idx_ok), max_slope(idx_ok,k), '-o', ...
         'Color', cmap2(k,:), 'LineWidth', 2, 'MarkerFaceColor', cmap2(k,:), ...
         'DisplayName', sprintf('%d rpm', rpm_values(k)));
end
set(gca, 'XScale', 'log');
xlabel('Fo  [m/s]', 'FontSize', 12);
ylabel('Max slope  [°]', 'FontSize', 12);
title('Maximum hh slope vs evaporation flux', 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 8, 'NumColumns', 2);
grid on;
hold off;
