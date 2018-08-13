function [componentTemplates,templateIDs] = ParseTemplates(sxStruct)
% This function parses a struct in SX format (converted from xml element)
% and returns a list of component templates in StructHA format.
% Additionally, a list of template names is returned, for fast name lookup.
% Order of templates reflects the xml file.
% References to other components are stored as indexes in the result array.

sxStruct = sxStruct.sspaceex{1};

% loop over all component definitions ("templates") in the xml
% process the contained data to convenient formats
% & store the resulting structs in "componentTemplates"
num_templates = length(sxStruct.component);
% preallocate cell array for templates
componentTemplates = cell(1,num_templates);

% store component names for fast resolution of references
templateIDs = strings(1,num_templates);

for i = 1:num_templates
    currentComp = sxStruct.component{i};
    % components can be "base" components, containing dynamics
    % they are identified by containing a "location"
    if isfield(currentComp,'location')
        fprintf("parsing BaseComponent %s...\n",currentComp.Attributes.id);
        parsedComp = ComputeBaseComponent(currentComp);
        parsedComp.isNetwork = false;
    % otherwise, components are "network" components, containing binds
    % binds reference other component templates
    elseif isfield(currentComp,'bind')
        fprintf("parsing NetworkComponent %s...\n",currentComp.Attributes.id);
        % represent references to other templates by their index in this array
        parsedComp = ComputeNetworkComponent(currentComp, templateIDs(1:i-1));
        parsedComp.isNetwork = true;
    else
        error("component #%i seems to be neither base nor network.",i);
    end
    
    componentTemplates{i} = parsedComp;
    templateIDs(i) = parsedComp.id;
end

fprintf("done parsing: %i components parsed\n",length(templateIDs));

end

