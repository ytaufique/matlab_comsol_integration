load('C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Results\rpm_sweep_spin_coat_1d.mat');
idx_ok = find(success);
rpm_ok = rpm_values(idx_ok);
hh_final_um = cellfun(@(v) v(end)*1e6, hh_cell(idx_ok));
figure; plot(rpm_ok, hh_final_um, 'k-^', 'LineWidth', 2, 'MarkerFaceColor', 'k');
xlabel('Rotation speed [rpm]'); ylabel('h at last radial point [µm]');
title('h value 50s vs rotation speed'); grid on;