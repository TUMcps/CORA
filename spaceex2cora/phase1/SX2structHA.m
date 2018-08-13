function [automaton,componentTemplates,componentInstances] = SX2structHA( xmlData , rootID )
%Input : 
%   xmlData: automaton decription file in SX format
%   name: Name of automaton object & filename of output matlab file
%   rootID: ID of SpaceEx component to be used as root component
%Output : automaton in structHA formate
%
%Example : SX2structHA('bball.xml','ball','IDball') returns bouncing ball
% structHA

disp("--------------------STEP 1 : PARSING XML FILE--------------------");
% compute Matlab structure of xml file
sxStruct = xml2struct(xmlData);

disp("--------------STEP 2 : PARSING COMPONENT DEFINITIONS-------------");
% parse component templates into individual structs
% structs are in StructHA format
[componentTemplates,templateIDs] = ParseTemplates(sxStruct);

% sanity check
if isempty(templateIDs)
    error("no component templates could be parsed");
end

% find index of "root" component
% store it in the field "rootIdx"
if nargin >= 2
    % if root ID is given, find the index of the corresponding template
    isRoot = string(rootID) == templateIDs;
    if ~any(isRoot)
        error("no components with the given rootID ""%s"" were found!",rootID);
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
% build the parallel automaton, by beginning at the root component
% and instantiating the tree of referenced components below
fprintf("building instance tree from root ""%s""...\n",rootID);
componentInstances = InstantiateComponents(componentTemplates,rootIdx);

disp("------------------STEP 4 : MERGING INSTANCE TREE-----------------");
% then merge the tree into a single BaseComponent
% (by computing the automaton product)
mergedComponent = ProductMerge(componentInstances);

disp("------------------STEP 5 : FORMALIZING AUTOMATON-----------------");
% quantize flow, invariant, guard & reset equations to matrices
formalizedComponent = FormalizeBaseComponent(mergedComponent);

% Package the automaton in the StructHA format
automaton.Components = {formalizedComponent};
% store name of source SX file 
[~,automaton.name,~] = fileparts(xmlData);
% store name of used root component
automaton.componentID = rootID;

disp("-----------------------SX2structHA COMPLETE-----------------------");

end

