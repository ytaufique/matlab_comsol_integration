% inspect_functions.m
% Extracts all analytic/interpolation function definitions from the model.

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try mphstart; catch; end

model = mphload('C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Project\spin_coat_stillwagon_topography_solvent_evaporation.mph');

outfile = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Matlab\spin_coat_evap_functions.txt';
fid = fopen(outfile, 'w');

props = {'funcname','expr','args','argunit','fununit', ...
         'filename','interp','extrap','smooth','pieces'};

ftags = model.func.tags();
fprintf(fid, 'FUNCTIONS (%d found)\n', length(ftags));
fprintf(fid, '%s\n', repmat('=',1,60));

for i = 1:length(ftags)
    tag   = char(ftags(i));
    label = char(model.func(tag).label());
    type  = class(model.func(tag));
    fprintf(fid, '\n[%s]  label: %s  type: %s\n', tag, label, type);
    for k = 1:length(props)
        p = props{k};
        try
            val = char(model.func(tag).getString(p));
            if ~isempty(strtrim(val))
                fprintf(fid, '  %-15s = %s\n', p, val);
            end
        catch; end
        try
            arr = model.func(tag).getStringArray(p);
            if ~isempty(arr)
                for m = 1:length(arr)
                    entry = char(arr(m));
                    if ~isempty(strtrim(entry))
                        fprintf(fid, '  %-12s[%d] = %s\n', p, m, entry);
                    end
                end
            end
        catch; end
    end
end

% Also dump dependent variable names from physics
fprintf(fid, '\n%s\nDEPENDENT VARIABLES\n%s\n', repmat('=',1,60), repmat('=',1,60));
try
    ptags = model.physics.tags();
    for i = 1:length(ptags)
        tag = char(ptags(i));
        try
            fields = model.physics(tag).field('dimensionless').component();
            fprintf(fid, '[%s] components: ', tag);
            for j = 1:length(fields)
                fprintf(fid, '%s  ', char(fields(j)));
            end
            fprintf(fid, '\n');
        catch e2
            fprintf(fid, '[%s] ERROR: %s\n', tag, e2.message);
        end
    end
catch e
    fprintf(fid, 'ERROR: %s\n', e.message);
end

fclose(fid);
fprintf('Done. Written to:\n  %s\n', outfile);
