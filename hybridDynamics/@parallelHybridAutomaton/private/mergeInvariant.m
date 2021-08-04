function res = mergeInvariant(obj, invList)
% mergeInvariant - Create the full dimensional invariant of the overall
%                  system from the invariants fo the subcomponents
%
% Syntax:  
%    res = mergeInvariant(obj, invList)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    invList - list with the invariant sets for all subcomponents
%
% Outputs:
%    res - resulting invariant set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Johann Schoepfer, Niklas Kochdumper
% Written:      08-June-2018  
% Last update:  09-July-2018 (NK, use "projectHighDim" function)
% Last revision: ---

%------------- BEGIN CODE --------------

    % project invariant set of the first subcomponent to high dimensional space 
    % of the overall automaton
    res = projectHighDim(invList{1},obj.numStates,obj.bindsStates{1});

    % loop over the invariants of the remaining subcomponents
    for i = 2:length(invList)

       % project set to high dimensional space of the overall automaton
       temp = projectHighDim(invList{i},obj.numStates,obj.bindsStates{i});

       % compute intersection with the invariants of the remaining
       % subcomponents
       res = res & temp;
    end
end

%------------- END OF CODE --------------