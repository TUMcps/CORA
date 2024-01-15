function res = test_zonoBundle_polytope
% test_zonoBundle_polytope - unit test function of polytope
%
% Syntax:
%    res = test_zonoBundle_polytope
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
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% fully-empty zonoBundle
% zB = zonoBundle.empty(2);
% P = polytope(zB);
% res(end+1,1) = representsa(P,'emptySet');

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
% convert to polytope
P = polytope(zB);
% vertices
V = vertices(zB);
V_ = vertices(P);
% compare results
res(end+1,1) = compareMatrices(V,V_,1e-12);

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
% convert to polytope
P = polytope(zB);
% vertices
V = vertices(zB);
V_ = vertices(P);
% compare results
res(end+1,1) = compareMatrices(V,V_,1e-12);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
