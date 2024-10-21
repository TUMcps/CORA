function invSet = mergeInvariants(pHA,locID,mergedLabels)
% mergeInvariants - creates the full-dimensional invariant of the overall
%    system from the invariants of the subcomponents
%
% Syntax:
%    invSet = mergeInvariants(pHA,locID,mergedLabels)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    locID - IDs of the currently active locations
%    mergedLabels - synchronization labels which are active in the
%                   composed location
%
% Outputs:
%    invSet - resulting invariant set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Johann Schoepfer, Niklas Kochdumper, Mark Wetzlinger
% Written:       08-June-2018  
% Last update:   09-July-2018 (NK, use "projectHighDim" function)
%                06-April-2023 (MW, integrate fullspace invariants)
%                24-May-2023 (MW, intersect invariants with transitions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize resulting invariant set
invSet = fullspace(pHA.nrOfStates);

% loop over the invariants of the remaining subcomponents
for i=1:length(pHA.components)

    % read out invariant of currently active location of i-th subcomponent
    inv = pHA.components(i).location(locID(i)).invariant;

    % read out set of outgoing transitions of currently active location
    transList = pHA.components(i).location(locID(i)).transition;

    % -> intersect invariant with complements of guard sets if the
    % corresponding transition's synchronization label is active in the
    % composed location
    for j=1:length(transList)
        if isa(transList(j).guard,'polytope') ...
                || isa(transList(j).guard,'levelSet') ...
                && ismember(transList(j).syncLabel,mergedLabels)
            % not-operation not implemented for some cases...
            try
                inv = and_(inv,~transList(j).guard,'exact');
            end
        end
    end

    % lift set to high-dimensional space of the overall automaton
    inv_lifted = lift_(inv,pHA.nrOfStates,pHA.bindsStates{i});

    % compute intersection with the invariants of the remaining
    % subcomponents
    invSet = and_(invSet,inv_lifted,'exact');

end

% ------------------------------ END OF CODE ------------------------------
