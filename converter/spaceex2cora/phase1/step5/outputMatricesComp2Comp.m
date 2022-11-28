function [C,D,k] = outputMatricesComp2Comp(states,inputs,outputs)
% outputMatricesComp2Comp - instantiates output matrices C, D, and output
%    vector k for output matrices of a component in a parallel hybrid
%    automaton; the outputs are inputs in other components, hence the
%    throughput matrix D and the constant offset k are all-zero (for now)
%
% Syntax:  
%    [C,D,k] = outputMatricesComp2Comp(states,inputs,outputs)
%
% Inputs:
%    states - all state variables
%    inputs - all input variables
%    outputs - all output variables
%
% Outputs:
%    C - output matrix
%    D - throughput matrix
%    k - output offset
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      25-June-2022
% Last update:  --- 
% Last revision:---

%------------- BEGIN CODE --------------

% if no outputs, then no equation y = Cx + Du + k
% hence, return empty matrices (different handling in data2parallelHA.m)
if isempty(outputs)
    C = []; D = []; k = [];
    return
end

% size init
sizeC = length(states);
sizeD = length(inputs);
sizek = length(outputs);
C = zeros(sizek,sizeC);
D = zeros(sizek,sizeD);
k = zeros(sizek,1);

% only C is potentially non-zero
for i=1:length(outputs)
    % see to which state the i-th output corresponds
    outputIdx = [states.name] == outputs(i).name;
    % add entry in C matrix
    C(i,outputIdx) = 1;
end

%------------- END OF CODE -------------