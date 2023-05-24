function res = test_parallelHybridAutomaton_isempty
% test_parallelHybridAutomaton_isempty - test function for isempty
%
% Syntax:  
%    res = test_parallelHybridAutomaton_isempty
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      19-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty parallel hybrid automata
res = isempty(parallelHybridAutomaton());

% full automaton
pHA = roomHeatingParallel();
res(end+1,1) = ~isempty(pHA);

% combine results
res = all(res);

%------------- END OF CODE --------------
