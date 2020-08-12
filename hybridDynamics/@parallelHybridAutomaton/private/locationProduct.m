function res = locationProduct(obj, loc)
% locationProduct - Construct a overall location object from the active 
%                   loctions of the subcomponts with a local automaton
%                   product
%
% Syntax:  
%    res = locationProduct(obj, loc)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    loc - ID's for the active location of each subcomponent
%
% Outputs:
%    res - constructed location object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Johann Schöpfer, Niklas Kochdumper
% Written:      08-June-2018  
% Last update:  09-July-2018
% Last revision: ---

%------------- BEGIN CODE --------------

    numComp = length(obj.components);

    % get active location objects
    locList = cell(1,numComp);

    for i = 1:numComp
        locList{i} = obj.components{i}.location{loc(i)};
    end

    % preallocate arrays
    invList = cell(1,numComp); 
    flowList = cell(1,numComp);
    transList = cell(1,numComp);

    % get properties from all active components
    for iLocation = 1:numComp
        invList{iLocation} = locList{iLocation}.invariant;
        flowList{iLocation} = locList{iLocation}.contDynamics;
        transList{iLocation} = locList{iLocation}.transition;
    end

    % merge invariants
    mergedInvSet = mergeInvariant(obj,invList);

    % merge flows (= continuous dynamics)
    mergedFlow = mergeFlows(obj,flowList);

    % merge transitions
    mergedTransSets = mergeTransitionSets(obj,transList,loc);

    % construct resulting location object
    res = location(mergedInvSet,mergedTransSets,mergedFlow);

end

%------------- END OF CODE --------------