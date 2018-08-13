function nc_out = ComputeNetworkComponent( nc_in , templateIDs )
%INPUT:
%   nc_in: component definition in SX format
%   (xml element converted to matlab struct)
%OUTPUT:
%   nc_out: network component in structHA format
%   contains:listOfVar (variable definitions, cont. states/inputs/params)
%            binds (included component templates (has idx,id & local name)
%               (also includes "mapping" for varialbe assignment)
%   CONVENTION: only strings leave this functions (no char arrays)

% Assign meta information.
nc_out.id = string(nc_in.Attributes.id);

% Call the function CollectVariables to CollectVariables and Constants.
% For internal use only, needed in function ComputeNetworkComponent and
% AssignTemplate
[listOfVar, ~] = CollectVariables(nc_in.param);
nc_out.listOfVar = listOfVar;
%nc_out.h_listOfLab = listOfLab; %unused since CORA doesn not support it

%preallocate data structure
num_binds = length(nc_in.bind);
nc_out.Binds = struct('idx',cell(1,num_binds));

for i = 1:num_binds
    
    % Find the id of the bound component template
    boundID = string(nc_in.bind{i}.Attributes.component);
    
    % search for this id in the list of already parsed templates
    % store this index for easier referencing of the template
    nc_out.Binds(i).idx = find(boundID == templateIDs,1);
    
    % store naming information
    nc_out.Binds(i).id = boundID;
    nc_out.Binds(i).localName = string(nc_in.bind{i}.Attributes.as);
    
    % parse variable mapping
    maps_in = nc_in.bind{i}.map;
    num_maps = length(maps_in);
    
    % preallocate arrays
    keys = strings(num_maps,1);
    values = sym(zeros(num_maps,1));
    values_text = strings(num_maps,1);
    for j = 1:num_maps
        % access name of the mapped variable & the mapped expression
        key_text = maps_in{j}.Attributes.key;
        value_text = maps_in{j}.Text;
        
        % Unfortunately, the symbolic Toolbox interprets variables named "i","j",
        % "I" or "J" as the imaginary number.
        % We perform a transformation on all variable names to avoid this.
        keys(j) = replaceImagVarnames(key_text);
        values_text(j) = replaceImagVarnames(value_text);
        
        % store mapped expression as a symbolic, to simplify variable
        % substitution in the future
        values(j) = str2symbolic(values_text(j));
    end
    nc_out.Binds(i).keys = keys;
    nc_out.Binds(i).values = values;
    nc_out.Binds(i).renames = values_text;
end


end


