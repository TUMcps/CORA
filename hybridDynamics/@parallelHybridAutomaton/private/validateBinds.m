function obj = validateBinds(obj)
% validateBinds - validity test for component interconnection
%
% This test ensures: -that inputBinds are valid
%
% Inputs:
%    obj - object of class parallelHybridAutomaton
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

    % construct stateBinds cell-array
    numComp = length(obj.components);
    stateBinds = cell(numComp,1);
    counter = 1;

    for i = 1:numComp
       numStates = obj.components{i}.location{1}.contDynamics.dim;
       stateBinds{i} = (counter:counter+numStates-1)';
       counter = counter + numStates;
    end
    
    obj.numStates = counter - 1;
    obj.bindsStates = stateBinds;

    % concatenate bind arrays vertically
    concatInputs = vertcat(obj.bindsInputs{:});

    % check if all indices are specified as integer values
    try
        int16(concatInputs);
    catch
        error('inputBinds: entries have to be integer values!');
    end

    % check whether inputBinds are valid
    inputComps = concatInputs(:,1);
    inputIndices = concatInputs(:,2);
    
    if any(inputComps < 0) || any(inputComps > numComp)
       error(['inputBinds: indices for components have to be greater ', ...
              'than 0 and small than the number of components!']); 
    end
    
    % get number of global inputs
    temp = unique(inputIndices(inputComps == 0));
    obj.numInputs = length(temp);

end

%------------- END OF CODE --------------

