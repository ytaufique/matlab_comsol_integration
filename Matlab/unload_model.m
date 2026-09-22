% unload_model.m
%
% Removes all models currently loaded in the COMSOL server from memory.

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try mphstart; catch; end
import com.comsol.model.util.*

tags = ModelUtil.tags();
if isempty(tags)
    fprintf('No models currently loaded.\n');
else
    for k = 1:length(tags)
        tag = char(tags(k));
        fprintf('Removing model: %s\n', tag);
        ModelUtil.remove(tag);
    end
    fprintf('Done. %d model(s) unloaded.\n', length(tags));
end
