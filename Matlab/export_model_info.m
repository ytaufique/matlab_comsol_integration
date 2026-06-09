% export_model_info.m
% Exports model parameters, variables, physics, mesh, materials, and
% study settings from a .mph file to a plain text file.

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try
    mphstart;
catch
end

disp('Loading model...');
model = mphload('C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Project\dip_coating_BUMP.mph');
disp('Model loaded.');

outfile = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Matlab\model_info.txt';
fid = fopen(outfile, 'w');
if fid == -1
    error('Could not open output file: %s', outfile);
end
disp('File opened for writing.');

writeline = @(s) fprintf(fid, '%s\n', s);
writesep  = @()  fprintf(fid, '%s\n', repmat('=', 1, 60));

writeline('COMSOL MODEL EXPORT');
writeline(datestr(now));
writesep();

disp('Writing parameters...');
%% Parameters
writeline('');
writeline('PARAMETERS');
writesep();
try
    names = model.param.varnames();
    for i = 1:length(names)
        name = char(names(i));
        val  = char(model.param.get(name));
        desc = char(model.param.descr(name));
        fprintf(fid, '  %-25s = %-25s  %s\n', name, val, desc);
    end
catch e
    writeline(['  ERROR: ' e.message]);
end

%% Variables
writeline('');
writeline('VARIABLES');
writesep();
try
    vtags = model.variable.tags();
    for i = 1:length(vtags)
        tag = char(vtags(i));
        fprintf(fid, '  Group: %s\n', tag);
        vnames = model.variable(tag).varnames();
        for j = 1:length(vnames)
            vname = char(vnames(j));
            vexpr = char(model.variable(tag).get(vname));
            vdesc = char(model.variable(tag).descr(vname));
            fprintf(fid, '    %-25s = %-30s  %s\n', vname, vexpr, vdesc);
        end
    end
catch e
    writeline(['  ERROR: ' e.message]);
end

%% Functions
writeline('');
writeline('FUNCTIONS');
writesep();
try
    ftags = model.func.tags();
    for i = 1:length(ftags)
        tag = char(ftags(i));
        label = char(model.func(tag).label());
        fprintf(fid, '  %-20s  label: %s\n', tag, label);
    end
catch e
    writeline(['  ERROR: ' e.message]);
end

%% Materials
writeline('');
writeline('MATERIALS');
writesep();
try
    mtags = model.material.tags();
    for i = 1:length(mtags)
        tag   = char(mtags(i));
        label = char(model.material(tag).label());
        fprintf(fid, '  %-20s  label: %s\n', tag, label);
    end
catch e
    writeline(['  ERROR: ' e.message]);
end

%% Geometry
writeline('');
writeline('GEOMETRY');
writesep();
try
    gtags = model.geom.tags();
    for i = 1:length(gtags)
        tag = char(gtags(i));
        fprintf(fid, '  tag: %s\n', tag);
    end
catch e
    writeline(['  ERROR: ' e.message]);
end

%% Mesh
writeline('');
writeline('MESH');
writesep();
try
    meshtags = model.mesh.tags();
    for i = 1:length(meshtags)
        tag = char(meshtags(i));
        fprintf(fid, '  tag: %s\n', tag);
        try
            stats = mphmeshstats(model, tag);
            fields = fieldnames(stats);
            for k = 1:length(fields)
                f = fields{k};
                v = stats.(f);
                if isnumeric(v)
                    fprintf(fid, '    %-25s = %s\n', f, mat2str(v));
                elseif ischar(v)
                    fprintf(fid, '    %-25s = %s\n', f, v);
                end
            end
        catch
        end
    end
catch e
    writeline(['  ERROR: ' e.message]);
end

%% Physics
writeline('');
writeline('PHYSICS');
writesep();

% Properties to attempt reading from each feature (covers General Form PDE
% and common BC types)
pde_props = {'f','Ga','ea','da','u','g','r','be','al','bt', ...
             'gamma','shape','order','units'};

try
    ptags = model.physics.tags();
    for i = 1:length(ptags)
        tag   = char(ptags(i));
        label = char(model.physics(tag).label());
        fprintf(fid, '\n  [%s]  %s\n', tag, label);
        try
            ftags = model.physics(tag).feature.tags();
            for j = 1:length(ftags)
                ftag   = char(ftags(j));
                flabel = char(model.physics(tag).feature(ftag).label());
                fprintf(fid, '    Feature: %-20s  (%s)\n', ftag, flabel);
                for k = 1:length(pde_props)
                    prop = pde_props{k};
                    try
                        val = model.physics(tag).feature(ftag).getString(prop);
                        val = char(val);
                        if ~isempty(val)
                            fprintf(fid, '      %-10s = %s\n', prop, val);
                        end
                    catch
                    end
                    try
                        arr = model.physics(tag).feature(ftag).getStringArray(prop);
                        if ~isempty(arr)
                            for m = 1:length(arr)
                                entry = char(arr(m));
                                if ~isempty(entry)
                                    fprintf(fid, '      %-10s[%d] = %s\n', prop, m, entry);
                                end
                            end
                        end
                    catch
                    end
                end
            end
        catch e2
            fprintf(fid, '    ERROR reading features: %s\n', e2.message);
        end
    end
catch e
    writeline(['  ERROR: ' e.message]);
end

%% Studies
writeline('');
writeline('STUDIES');
writesep();

study_props = {'tlist', 'tunit', 'tols', 'rtol', 'atol', 'dtmax', 'dtmin', 'timestepmax'};

try
    stags = model.study.tags();
    for i = 1:length(stags)
        tag   = char(stags(i));
        label = char(model.study(tag).label());
        fprintf(fid, '\n  [%s]  %s\n', tag, label);
        try
            ftags = model.study(tag).feature.tags();
            for j = 1:length(ftags)
                ftag   = char(ftags(j));
                flabel = char(model.study(tag).feature(ftag).label());
                fprintf(fid, '    Feature: %-20s  (%s)\n', ftag, flabel);
                for k = 1:length(study_props)
                    prop = study_props{k};
                    try
                        val = char(model.study(tag).feature(ftag).getString(prop));
                        if ~isempty(val)
                            fprintf(fid, '      %-15s = %s\n', prop, val);
                        end
                    catch
                    end
                end
            end
        catch e2
            fprintf(fid, '    ERROR: %s\n', e2.message);
        end
    end
catch e
    writeline(['  ERROR: ' e.message]);
end

fclose(fid);
fprintf('Done. Exported to:\n  %s\n', outfile);
