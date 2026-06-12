% inspect_spin_coat_evap.m
% Exports full model info from spin_coat_stillwagon_topography_solvent_evaporation.mph
% Run this once so we can read the PDE formulation before building the 2D version.

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try
    mphstart;
catch
end

model_path = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Project\spin_coat_stillwagon_topography_solvent_evaporation.mph';
disp('Loading model...');
model = mphload(model_path);
disp('Model loaded.');

outfile = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Matlab\spin_coat_evap_info.txt';
fid = fopen(outfile, 'w');

writeline = @(s) fprintf(fid, '%s\n', s);
writesep  = @()  fprintf(fid, '%s\n', repmat('=', 1, 60));

writeline('spin_coat_stillwagon_topography_solvent_evaporation.mph');
writeline(datestr(now));
writesep();

%% Parameters
writeline(''); writeline('PARAMETERS'); writesep();
try
    names = model.param.varnames();
    for i = 1:length(names)
        name = char(names(i));
        val  = char(model.param.get(name));
        desc = char(model.param.descr(name));
        fprintf(fid, '  %-25s = %-30s  %s\n', name, val, desc);
    end
catch e; writeline(['  ERROR: ' e.message]); end

%% Variables
writeline(''); writeline('VARIABLES'); writesep();
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
            fprintf(fid, '    %-25s = %-35s  %s\n', vname, vexpr, vdesc);
        end
    end
catch e; writeline(['  ERROR: ' e.message]); end

%% Geometry
writeline(''); writeline('GEOMETRY'); writesep();
try
    gtags = model.geom.tags();
    for i = 1:length(gtags)
        tag = char(gtags(i));
        sdim = model.geom(tag).getSDim();
        fprintf(fid, '  tag: %-10s  sdim: %d\n', tag, sdim);
    end
catch e; writeline(['  ERROR: ' e.message]); end

%% Mesh
writeline(''); writeline('MESH'); writesep();
try
    meshtags = model.mesh.tags();
    for i = 1:length(meshtags)
        tag = char(meshtags(i));
        try
            stats = mphmeshstats(model, tag);
            fprintf(fid, '  tag: %s  sdim=%d  nelems=%s\n', tag, stats.sdim, mat2str(stats.numelem));
        catch
            fprintf(fid, '  tag: %s\n', tag);
        end
    end
catch e; writeline(['  ERROR: ' e.message]); end

%% Physics — full property dump
writeline(''); writeline('PHYSICS'); writesep();
all_props = {'f','Ga','ea','da','u','g','r','be','al','bt','gamma', ...
             'shape','order','units','evap','E','E_evap','ke','c','flux'};
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
                fprintf(fid, '    Feature: %-25s  (%s)\n', ftag, flabel);
                for k = 1:length(all_props)
                    prop = all_props{k};
                    try
                        arr = model.physics(tag).feature(ftag).getStringArray(prop);
                        for m = 1:length(arr)
                            entry = char(arr(m));
                            if ~isempty(strtrim(entry))
                                fprintf(fid, '      %-10s[%d] = %s\n', prop, m, entry);
                            end
                        end
                    catch; end
                    try
                        val = char(model.physics(tag).feature(ftag).getString(prop));
                        if ~isempty(strtrim(val))
                            fprintf(fid, '      %-10s    = %s\n', prop, val);
                        end
                    catch; end
                end
            end
        catch e2
            fprintf(fid, '    ERROR: %s\n', e2.message);
        end
    end
catch e; writeline(['  ERROR: ' e.message]); end

%% Study
writeline(''); writeline('STUDIES'); writesep();
study_props = {'tlist','tunit','rtol','atol','dtmax','dtmin'};
try
    stags = model.study.tags();
    for i = 1:length(stags)
        tag   = char(stags(i));
        label = char(model.study(tag).label());
        fprintf(fid, '\n  [%s]  %s\n', tag, label);
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
                catch; end
            end
        end
    end
catch e; writeline(['  ERROR: ' e.message]); end

fclose(fid);
fprintf('Done. Exported to:\n  %s\n', outfile);
