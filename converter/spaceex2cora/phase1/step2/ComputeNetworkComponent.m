function nc_out = ComputeNetworkComponent(nc_in,templateIDs)
% ComputeNetworkComponent - For a network component of the SpaceEx model,
%    extract properties of the component (mapping of variables, connection
%    between components etc.) from the struct and prepare them for later
%    conversion to CORA format
%    note: only strings leave this function (no char arrays)!
%
% Syntax:
%    nc_out = ComputeNetworkComponent(nc_in,templateIDs)
%
% Inputs:
%    nc_in (struct) - component definition in SX format
%                     (xml element converted to matlab struct)
%    templateIDs - list of component names
%
% Outputs:
%    nc_out (struct) - network component in structHA format, with fields
%            .listOfVar (variable definitions, cont. states/inputs/params)
%            .binds (included component templates (has idx,id & local name,
%                    also includes "mapping" for variable assignment)
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Assign meta information
nc_out.id = string(nc_in.Attributes.id);

% Call the function CollectVariables to collect variables and constants.
[listOfVars, listOfLabels] = CollectVariables(nc_in.param);
nc_out.listOfVar = listOfVars;
nc_out.listOfLabels = listOfLabels;

% init data structure
num_binds = length(nc_in.bind);
nc_out.Binds = struct('idx',cell(1,num_binds));

% loop over all binds (= base components that are instantiated in the given
% network component)
for i = 1:num_binds
    
    % Find the id (= name) of the bound base component template
    boundID = string(nc_in.bind{i}.Attributes.component);
    
    % search for this id in the list of already parsed templates and store
    % this index for easier referencing of the template
    nc_out.Binds(i).idx = find(boundID == templateIDs,1);
    
    % store naming information
    nc_out.Binds(i).id = boundID;
    nc_out.Binds(i).localName = string(nc_in.bind{i}.Attributes.as);
    
    
    % parse variable mapping
    maps_in = nc_in.bind{i}.map;
    
    % delete synchronization labels from mapping before bind-generation
    % (sychronization labels here lead to errors due to indexing later)
    labelList = arrayfun(@(x) x.name,listOfLabels,'UniformOutput',false);
    keyList = cellfun(@(x) x.Attributes.key,maps_in,'UniformOutput',false);
    maps_in = maps_in(~ismember(keyList,labelList));
    
    % pre-allocate arrays
    num_maps = length(maps_in);
    keys = strings(num_maps,1);
    values = sym(zeros(num_maps,1));
    values_text = strings(num_maps,1);
    
    for j = 1:num_maps
        % name of the mapped variable & the mapped expression
        key_text = maps_in{j}.Attributes.key;
        value_text = maps_in{j}.Text;
        
        % Unfortunately, the Symbolic Toolbox interprets variables named
        % "i","j", "I" or "J" as the imaginary number.
        % -> transform all variable names to avoid this.
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

% ------------------------------ END OF CODE ------------------------------
