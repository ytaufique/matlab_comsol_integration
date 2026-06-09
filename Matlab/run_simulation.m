% run_simulation.m
% Runs the dip coating simulation and plots hh and pp at t_end.
% Supports optional parameter sweep over one or more model parameters.

% -------------------------------------------------------------------------
% TIME SETTINGS
% use_model_time = true  : use time range from .mph file
% use_model_time = false : override with t_start, t_end, dt below
use_model_time = true;
t_start = 0;       % seconds
t_end   = 1;       % seconds
dt      = 0.01;    % seconds
% -------------------------------------------------------------------------
% PARAMETER SWEEP
% do_sweep = false : single run with model/params.json values
% do_sweep = true  : loop over all sweep entries below
%
% sweep_mode:
%   'grid'   - run all combinations of all parameter values
%   'paired' - vary all parameters together (arrays must be same length)
%
% Add or remove sweeps(n) blocks as needed.
% Check params.json for exact parameter names.
do_sweep   = true;
sweep_mode = 'grid';

sweeps(1).param  = 'mu';
sweeps(1).values = [0.001, 0.005, 0.01];
sweeps(1).unit   = '[Pa*s]';

sweeps(2).param  = 'U';
sweeps(2).values = [0.003, 0.005];
sweeps(2).unit   = '[m/s]';
% -------------------------------------------------------------------------

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try
    mphstart;
catch
end

%% Load model
disp('Loading model...');
model = mphload('C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Project\dip_coating_BUMP_1.mph');

%% Apply baseline parameters from params.json if it exists
paramsfile = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Matlab\params.json';
if isfile(paramsfile)
    disp('Loading parameters from params.json...');
    params = jsondecode(fileread(paramsfile));
    pfields = fieldnames(params);
    for i = 1:length(pfields)
        model.param.set(pfields{i}, params.(pfields{i}).value);
    end
end

%% Set time list
if use_model_time
    tlist = char(model.study('std1').feature('time').getString('tlist'));
    fprintf('Using time list from model: %s\n', tlist);
else
    tlist = sprintf('range(%g,%g,%g)', t_start, dt, t_end);
    model.study('std1').feature('time').set('tlist', tlist);
    fprintf('Using custom time list: %s\n', tlist);
end

%% Build list of parameter combinations to run
if do_sweep
    n_params = length(sweeps);
    if strcmp(sweep_mode, 'grid')
        % Build index grid for all combinations
        val_lengths = arrayfun(@(s) length(s.values), sweeps);
        idx_arrays  = arrayfun(@(n) 1:n, val_lengths, 'UniformOutput', false);
        grids       = cell(1, n_params);
        [grids{:}]  = ndgrid(idx_arrays{:});
        n_runs      = numel(grids{1});
        run_indices = zeros(n_runs, n_params);
        for p = 1:n_params
            run_indices(:, p) = grids{p}(:);
        end
    else % paired
        n_runs      = length(sweeps(1).values);
        run_indices = repmat((1:n_runs)', 1, n_params);
    end
else
    n_runs      = 1;
    run_indices = [];
end

%% Pre-create sweep figures
fig_hh = figure; hold on;
xlabel('x (mm)'); ylabel('Film thickness h (\mum)');
grid on;

fig_pp = figure; hold on;
xlabel('x (mm)'); ylabel('Pressure p (Pa)');
grid on;

if do_sweep
    title(fig_hh.CurrentAxes, sprintf('Film Thickness at t_{end} (%s sweep)', sweep_mode));
    title(fig_pp.CurrentAxes, sprintf('Pressure at t_{end} (%s sweep)', sweep_mode));
else
    title(fig_hh.CurrentAxes, 'Film Thickness at t_{end}');
    title(fig_pp.CurrentAxes, 'Pressure at t_{end}');
end

sweep_colors = lines(n_runs);

import com.comsol.model.util.ModelUtil;
ModelUtil.showProgress(true);

%% Run loop
for k = 1:n_runs

    % Set parameters for this run
    label_parts = cell(1, size(run_indices, 2));
    if do_sweep
        for p = 1:n_params
            idx    = run_indices(k, p);
            sv     = sweeps(p).values(idx);
            sval   = sprintf('%g%s', sv, sweeps(p).unit);
            model.param.set(sweeps(p).param, sval);
            label_parts{p} = sprintf('%s=%g', sweeps(p).param, sv);
        end
        run_label = strjoin(label_parts, ', ');
        fprintf('\n--- Run %d/%d: %s ---\n', k, n_runs, run_label);
    else
        run_label = 'default';
    end

    % Run study
    model.study('std1').run();
    fprintf('Run %d complete.\n', k);

    % Plot COMSOL result groups on first run only
    if k == 1
        pgtags = model.result.tags();
        for i = 1:length(pgtags)
            try
                mphplot(model, char(pgtags(i)), 'rangenum', 1);
            catch
            end
        end
    end

    % Extract hh and pp at final time
    solinfo = mphsolinfo(model, 'soltag', 'sol1');
    idx_end = length(solinfo.solvals);
    d_hh    = mpheval(model, 'hh', 'edim', 1, 'solnum', idx_end);
    d_pp    = mpheval(model, 'pp', 'edim', 1, 'solnum', idx_end);
    xcoord  = d_hh.p(1,:);

    % Add to sweep figures
    figure(fig_hh);
    plot(xcoord * 1e3, d_hh.d1 * 1e6, 'Color', sweep_colors(k,:), ...
        'LineWidth', 1.5, 'DisplayName', run_label);

    figure(fig_pp);
    plot(xcoord * 1e3, d_pp.d1, 'Color', sweep_colors(k,:), ...
        'LineWidth', 1.5, 'DisplayName', run_label);

end

figure(fig_hh); legend('Location', 'best');
figure(fig_pp); legend('Location', 'best');
