function res = test_zonotope_isemptyobject
% test_zonotope_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_zonotope_isemptyobject
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
% Written:       03-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate zonotopes
Z1 = zonotope.empty(2);
Z2 = zonotope([1;1]);
Z3 = zonotope([1;1],[1 3 -2; 2 -4 2]);

% check results
res = isemptyobject(Z1) && ~isemptyobject(Z2) && ~isemptyobject(Z3);

% ------------------------------ END OF CODE ------------------------------
