% inspect_model.m
% Loads the dip coating model and prints all tags and parameters.
% Run this first to understand the model structure before building
% a full simulation script.

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try
    mphstart;
catch
end

model = mphload('C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Project\dip_coating_BUMP.mph');

disp('--- Study tags ---')
disp(model.study.tags())

disp('--- Physics tags ---')
disp(model.physics.tags())

disp('--- Mesh tags ---')
disp(model.mesh.tags())

disp('--- Parameters ---')
names = model.param.varnames();
for i = 1:length(names)
    name = char(names(i));
    val  = char(model.param.get(name));
    desc = char(model.param.descr(name));
    fprintf('  %-20s = %-20s  (%s)\n', name, val, desc);
end
