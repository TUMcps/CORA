function pHA = validateBinds(pHA)
% validateBinds - validity test for component interconnection (inputBinds)
%
% Syntax:
%    pHA - validateBinds(pHA)
% 
% Inputs:
%    pHA - parallelHybridAutomaton object
%    inputBinds - cell array of mx2 int arrays. Maps component inputs to 
%                 states/inputs of composed system
%
% Outputs:
%    pHA - parallelHybridAutomaton object (properties 'bindsStates',
%          'numInputs', and 'numStates' now set)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Johann Schoepfer
% Written:      05-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% construct stateBinds cell-array
numComp = length(pHA.components);
stateBinds = cell(numComp,1);
counter = 1;

for i=1:numComp
    numStates = pHA.components{i}.location{1}.contDynamics.dim;
    stateBinds{i} = (counter:counter+numStates-1)';
    counter = counter + numStates;
end

pHA.numStates = counter - 1;
pHA.bindsStates = stateBinds;

% concatenate bind arrays vertically
concatInputs = vertcat(pHA.bindsInputs{:});

% check whether inputBinds are valid
inputComps = concatInputs(:,1);
inputIndices = concatInputs(:,2);

% get number of global inputs
temp = unique(inputIndices(inputComps == 0));
pHA.numInputs = length(temp);

%------------- END OF CODE --------------
