function res = test_simResult_isemptyobject
% test_simResult_isemptyobject - unit test function for isemptyobject
%
% Syntax:
%    res = test_simResult_isemptyobject()
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
% Written:       01-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty simResult
simRes = simResult();
assert(isemptyobject(simRes));

% simResult with trajectory
t = {[0; 0.02; 0.05]};
x = {[1 1; 0.9 1.1; 0.8 1.2]};
simRes = simResult(x,t);
assert(~isemptyobject(simRes));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
