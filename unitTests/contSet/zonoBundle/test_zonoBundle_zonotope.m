function res = test_zonoBundle_zonotope
% test_zonoBundle_zonotope - unit test function of zonotope conversion
%
% Syntax:
%    res = test_zonoBundle_zonotope
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

% fully-empty zonoBundle
zB = zonoBundle.empty(2);
Z = zonotope(zB);
res = representsa(Z,'emptySet');

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
% convert to zonotope
Z = zonotope(zB);
res(end+1,1) = contains(Z,zB);

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
% convert to zonotope
Z = zonotope(zB);
res(end+1,1) = contains(Z,zB);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
