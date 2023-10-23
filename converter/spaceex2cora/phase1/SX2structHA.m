function [automaton,componentTemplates,componentInstances] = ...
    SX2structHA(xmlData,convtype,rootID,name)
% SX2structHA - reads a SpaceEx xml-file and returns a struct file of a
%     (parallel or flat) hybrid automaton
%
% Syntax:
%    [automaton,componentTemplates,componentInstances] = ...
%       SX2structHA(xmlData,convtype)
%    [automaton,componentTemplates,componentInstances] = ...
%       SX2structHA(xmlData,convtype,rootID)
%    [automaton,componentTemplates,componentInstances] = ...
%       SX2structHA(xmlData,convtype,rootID,name)
%
% Inputs:
%    xmlData - automaton description file in SpaceEx format (xml-file)
%    convtype - type of conversion (parallel HA = true, flat HA = false)
%    rootID - (optional) ID of SpaceEx component which should be used as
%             the root component
%    name - (optional) name of automaton object & filename of output matlab
%           file
%
% Outputs:
%    automaton (struct) - all information about the hybrid automaton
%                         .components: list of components
%                         .name: name of automaton (same as input arg)
%                         .componentID: id number of root component
%    componentTemplates - list of all base/network components
%    componentInstances - processed list of all base/network components
%
% Example: 
%    % requires file 'bball.xml' in cora/models/SpaceEx
%    SX2structHA('bball.xml',false,'IDball','ball');
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

if convtype
    disp("--------------------CONVERSION TO PARALLEL HA--------------------");
else
    disp("----------------------CONVERSION TO FLAT HA----------------------");
end

disp("--------------------STEP 1 : PARSING XML FILE--------------------");
% compute Matlab structure of xml file
sxStruct = xml2struct(xmlData);
% sxStruct.sspaceex is a struct with fields
% - Attributes (meta information about SpaceEx)
% - component (cell-array of base/network components = BC/NC)
%   BC have the following fields:
%   - param (cell-array of structs): with fields
%     .Text = ""
%     .Attributes (struct) = Information about variable, i.e., name,
%     controlled, local, dynamics, type, ...
%   - location (cell-array of structs): with fields
%     .invariant (cell-array of structs): .Text = invariant
%     .flow (cell-array of structs): .Text = flow equation
%     .Attributes (struct): meta-data about location (position, name, ...)
%   - transition (cell-array of structs): with fields
%     .label (cell-array of structs): .Text = synchronization label
%     .assignment (cell-array of structs): .Text = reset function
%     .labelposition: (cell-array of structs): meta-data about label
%     .middlepoint (cell-array of structs): meta-data about transition
%     .Attributes (struct): .source/.target = indices
%   - Attributes (struct): .id = name of BC
%   NC have the following fields:
%   - param (cell-array of structs): .Attributes = Information about
%     variable, i.e., name, controlled, local, dynamics, type, ...
%   - bind (cell-array of structs): with fields
%     .map (cell-array of structs): with fields
%      .Text = name of variable
%      .Attributes (struct): .key = bind of variable
%     .Attributes (struct): meta-data about network component
%   - Attributes (struct): .id = name of NC


disp("--------------STEP 2 : PARSING COMPONENT DEFINITIONS-------------");
% parse component templates into individual structs
[componentTemplates,templateIDs] = ParseTemplates(sxStruct);
% componentTemplates contains all BC and NC
% templateIDs contains the names of all BC/NC (for quick lookup)

% sanity check: there has to be at least one parsed component
if isempty(templateIDs)
    throw(CORAerror('CORA:converterIssue',...
        'No component templates could be parsed'));
end

% find the index (rootIdx) and name (rootID) of the "root" component
if nargin >= 3 && ~isempty(rootID)
    % if rootID is given, find the index of the corresponding template
    isRoot = string(rootID) == templateIDs;
    if ~any(isRoot)
        % given rootID not found in any BC/NC
        throw(CORAerror('CORA:converterIssue',...
            ['No components with the given rootID ' char(rootID) ' were found!']));
    else
        % use index of first search hit
        rootIdx = find(isRoot,1);
    end
else
    % no root ID is given, default last item as root
    rootIdx = length(componentTemplates);
    rootID = componentTemplates{rootIdx}.id;
end


disp("--------------STEP 3 : BUILDING COMPONENT INSTANCES--------------");
% build the (parallel) automaton, by beginning at the root component and
% instantiating the tree of 'children' components below
fprintf("building instance tree from root '%s'...\n",rootID);
componentInstances = InstantiateComponents(componentTemplates,rootIdx);


if ~convtype
    % flat hybrid automaton
    disp("------------------STEP 4 : MERGING INSTANCE TREE-----------------");
    % merge the tree into a single base component by computing the
    % automaton product
    mergedComponent = ProductMerge(componentInstances);
else
    % parallel hybrid automaton
    disp("-----------------STEP 4 : BUILDING COMPONENT LIST----------------");
    % manipulate base components such that the function
    % FormalizeBaseComponent (see step 5) can be used
    componentList = createBClist(componentInstances);
end


disp("------------------STEP 5 : FORMALIZING AUTOMATON-----------------");
% rewrite symbolic flow, invariant, guard, and reset equations as matrices
if ~convtype
    % determine states, inputs, and outputs from merged base component
    mergedComponent = classifyVariablesFlat(mergedComponent);
    % flat hybrid automaton: only one component to be formalized
    % package the automaton in the StructHA format
    automaton.components = {FormalizeBaseComponent(mergedComponent,convtype)};
else    
    % parallel hybrid automaton
    % go over all base components to check which variables are states,
    % inputs, and outputs in which component
    componentList = classifyVariablesParallel(componentList);
    % loop over parallel (base) components
    for i=1:length(componentList)
        % package the automaton in the StructHA format
        automaton.components(i) = {FormalizeBaseComponent(componentList(i),convtype)};
    end
end

% store name of source SpaceEx file 
if nargin == 4
    automaton.name = name;
else
    [~,automaton.name,~] = fileparts(xmlData);
end
% store name of used root component
automaton.componentID = rootID;

disp("-----------------------SX2structHA COMPLETE----------------------");

end

% ------------------------------ END OF CODE ------------------------------
