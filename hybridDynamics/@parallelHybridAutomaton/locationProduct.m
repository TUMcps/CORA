function loc = locationProduct(pHA,locID,allLabels)
% locationProduct - Constructs an overall location object from the active 
%    loctions of the subcomponents with a local automaton product
%
% Syntax:
%    loc = locationProduct(pHA,locID,allLabels)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    locID - IDs of the currently active locations
%    allLabels - information about synchronization labels
%
% Outputs:
%    loc - constructed location object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Johann Schoepfer, Niklas Kochdumper
% Written:       08-June-2018  
% Last update:   09-July-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of components in parallel hybrid automaton
numComp = length(pHA.components);

% preallocate arrays for merging of invariants, flows, and transition sets
invList = cell(1,numComp); 
flowList = cell(1,numComp);
transList = cell(1,numComp);

for i=1:numComp
    % read out active location of i-th subcomponent
    curr_loc = pHA.components(i).location(locID(i));
    % write invariant, flow, transition set into array for merge
    invList{i} = curr_loc.invariant;
    flowList{i} = curr_loc.contDynamics;
    transList{i} = curr_loc.transition;
end

% merge transition sets
[mergedTransSets,mergedLabels] = mergeTransitionSets(pHA,transList,locID,allLabels);

% merge invariants
mergedInvSet = mergeInvariants(pHA,invList,transList,mergedLabels);

% merge flows
mergedFlow = mergeFlows(pHA,flowList,locID);

% construct resulting location object
loc = location(mergedInvSet,mergedTransSets,mergedFlow);

% ------------------------------ END OF CODE ------------------------------
