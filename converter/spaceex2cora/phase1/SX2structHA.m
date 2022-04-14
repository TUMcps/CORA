function [automaton,componentTemplates,componentInstances] = SX2structHA(xmlData, convtype, rootID, name)
% Input : 
%   xmlData: automaton description file in SX format
%   convtype: type of conversion (parallel HA = 1/ flat HA = 0)
%   rootID: ID of SpaceEx component to be used as root component
%   name: Name of automaton object & filename of output matlab file
% Output : automaton in structHA format
%
% Example : SX2structHA('bball.xml',0,'IDball','ball')
%   returns structHA (will become flat HA) of bball.xml

%------------- BEGIN CODE --------------

if convtype
    disp("--------------------CONVERSION TO PARALLEL HA--------------------");
else
    disp("----------------------CONVERSION TO FLAT HA----------------------");
end

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
if nargin >= 3 && ~isempty(rootID)
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

if ~convtype
    disp("------------------STEP 4 : MERGING INSTANCE TREE-----------------");
    % FLAT HA:
    % then merge the tree into a single BaseComponent
    % (by computing the automaton product)
    mergedComponent = ProductMerge(componentInstances);
else
    disp("-----------------STEP 4 : BUILDING COMPONENT LIST----------------");
    % PARALLEL HA:
    % manipulate BaseComponents such that FormalizeBaseComponent can be used
    componentList = createBClist(componentInstances);
end
    
disp("------------------STEP 5 : FORMALIZING AUTOMATON-----------------");
% quantize flow, invariant, guard & reset equations to matrices
if ~convtype
    formalizedComponent = FormalizeBaseComponent(mergedComponent,convtype);
    % Package the automaton in the StructHA format
    automaton.Components = {formalizedComponent};
else    
    % loop over parallel components
    for i=1:size(componentList,2)
        formalizedComponents(i) = FormalizeBaseComponent(componentList(i),convtype);
        % Package the automaton in the StructHA format
        automaton.Components(i) = {formalizedComponents(i)};
    end
end

% store name of source SX file 
if nargin == 4
    automaton.name = name;
else
    [~,automaton.name,~] = fileparts(xmlData);
end
% store name of used root component
automaton.componentID = rootID;

disp("-----------------------SX2structHA COMPLETE----------------------");

end

%------------- END OF CODE --------------