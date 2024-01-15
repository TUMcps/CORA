function res = test_zonoBundle_dim
% test_zonoBundle_dim - unit test function of dim
%
% Syntax:
%    res = test_zonoBundle_dim
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
n = 2;
zB = zonoBundle.empty(n);
res = dim(zB) == 2;

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
n = dim(zB);
res(end+1,1) = n == 2;

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
n = dim(zB);
res(end+1,1) = n == 2;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
