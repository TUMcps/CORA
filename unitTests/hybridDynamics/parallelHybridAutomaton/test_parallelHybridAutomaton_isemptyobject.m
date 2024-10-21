function res = test_parallelHybridAutomaton_isemptyobject
% test_parallelHybridAutomaton_isemptyobject - test function for
%    isemptyobject
%
% Syntax:
%    res = test_parallelHybridAutomaton_isemptyobject
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

% Authors:       Mark Wetzlinger
% Written:       19-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty parallel hybrid automata
assert(isemptyobject(parallelHybridAutomaton()));

% full automaton
pHA = roomHeatingParallel();
assert(~isemptyobject(pHA));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
