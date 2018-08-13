function [numStates,numInputs] = validateBinds(stateBinds,inputBinds)
% validateBinds - validity test for component interconnection
%
% This test ensures: -that component stateBinds are a partition
%                    -that inputBinds are valid
% Inputs:
%    stateBinds - cell array of nx1 integer arrays. Maps component states 
%                 to states of composed system
%    inputBinds - cell array of nx2 int arrays. Maps component inputs to 
%                 states/inputs of composed system
%
% Outputs:
%    numStates - number of states of the composed system
%    numInputs - number of inputs of the composed system
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Johann Schoepfer
% Written: 05-June-2018
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % concatenate bind arrays vertically
    concatStates = vertcat(stateBinds{:});
    concatInputs = vertcat(inputBinds{:});

    % check if all indices are specified as integer values
    try
        int16(concatStates);
        int16(concatInputs);
    catch
        error('stateBinds/inputBinds: entries have to be integer values!');
    end

    % check whether the stateBinds form a partition of continuous indices
    uniqueStates = unique(concatStates);
    numStates = length(uniqueStates);

    if length(concatStates) ~= numStates
        error('stateBinds: two states were bound to the same variable!');
    elseif any(uniqueStates<1) || any(uniqueStates>numStates)
        error('stateBinds: must map to a continuous index beginning at 1!');
    end

    % check whether inputBinds are valid
    inputComps = concatInputs(:,1);
    inputIndices = concatInputs(:,2);
    
    if any(inputComps < 0)
       error('inputBinds: indices for components have to be greater than 0!'); 
    end
    
    % get number of global inputs
    ind = find(inputComps == 0);
    temp = unique(inputIndices(ind));
    numInputs = length(temp);

end

%------------- END OF CODE --------------

