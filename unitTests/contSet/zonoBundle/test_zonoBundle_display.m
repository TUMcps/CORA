function res = test_zonoBundle_display
% test_zonoBundle_display - unit test function of display
%
% Syntax:
%    res = test_zonoBundle_display
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% fully-empty zonoBundle
zB = zonoBundle.empty(2)

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2})

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2})

% empty generator matrix
Z3 = zonotope([0;1]);
zB = zonoBundle({Z1,Z3})

% ------------------------------ END OF CODE ------------------------------
