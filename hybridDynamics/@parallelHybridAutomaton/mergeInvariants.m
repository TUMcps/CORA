function invSet = mergeInvariants(pHA,invList)
% mergeInvariants - creates the full-dimensional invariant of the overall
%    system from the invariants of the subcomponents; all subcomponents
%    with an empty invariant (= full space) are skipped
%
% Syntax:  
%    invSet = mergeInvariants(pHA,invList)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    invList - list with the invariant sets for all subcomponents
%
% Outputs:
%    invSet - resulting invariant set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Johann Schoepfer, Niklas Kochdumper
% Written:      08-June-2018  
% Last update:  09-July-2018 (NK, use "projectHighDim" function)
% Last revision:---

%------------- BEGIN CODE --------------

% initialize resulting invariant set
invSet = [];

% loop over the invariants of the remaining subcomponents
for i=1:length(invList)

    % project set to high dimensional space of the overall automaton
    if ~(isnumeric(invList{i}) && isempty(invList{i}))
        temp = projectHighDim(invList{i},pHA.numStates,pHA.bindsStates{i});

        % compute intersection with the invariants of the remaining
        % subcomponents
        % note: if an invariant is empty, we assume it stems from a
        % converted model which had no invariant set (empty invariant)
        % in the current location; this poses no problem for further
        % analysis since it is only checked whether the reachable set fully
        % exits the invariant set, for which the remaining full-dimensional
        % invariants of the other subcomponents suffice
        if isnumeric(invSet) && isempty(invSet)
            % first full-dimensional invariant
            invSet = temp;
        else
            invSet = invSet & temp;
        end

    end

end

%------------- END OF CODE --------------