function invSet = mergeInvariants(pHA,invList,transList,mergedLabels)
% mergeInvariants - creates the full-dimensional invariant of the overall
%    system from the invariants of the subcomponents
%
% Syntax:
%    invSet = mergeInvariants(pHA,invList,transList,mergedLabels)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    invList - list with the invariant sets for all subcomponents
%    transList - list with the transition for all subcomponents
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
invSet = fullspace(pHA.dim);

% loop over the invariants of the remaining subcomponents
for i=1:length(invList)

    % intersect invariant with complements of guard sets if the
    % corresponding transition's synchronization label is active in the
    % composed location
    for j=1:length(transList{i})
        if isa(transList{i}(j).guard,'polytope') ...
                || isa(transList{i}(j).guard,'levelSet')...
                && ismember(transList{i}(j).syncLabel,mergedLabels)
            % not-operation not implemented for some cases...
            try
                invList{i} = and_(invList{i},~transList{i}(j).guard,'exact');
            end
        end
    end

    % lift set to high-dimensional space of the overall automaton
    temp = lift_(invList{i},pHA.dim,pHA.bindsStates{i});

    % compute intersection with the invariants of the remaining
    % subcomponents
    invSet = and_(invSet,temp,'exact');

end

% ------------------------------ END OF CODE ------------------------------
