% export_params.m
% Exports all model parameters to params.json for editing in VS Code.
% Run this once to generate the file, then edit params.json and run
% run_simulation.m to apply the changes.

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try
    mphstart;
catch
end

model = mphload('C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Project\dip_coating_BUMP.mph');

names = model.param.varnames();
params = struct();
for i = 1:length(names)
    name = char(names(i));
    val  = char(model.param.get(name));
    desc = char(model.param.descr(name));
    params.(name).value       = val;
    params.(name).description = desc;
end

outfile = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Matlab\params.json';
fid = fopen(outfile, 'w');
fprintf(fid, '%s', jsonencode(params, 'PrettyPrint', true));
fclose(fid);

fprintf('Parameters exported to:\n  %s\n', outfile);
