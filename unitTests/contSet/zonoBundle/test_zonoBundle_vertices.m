function res = test_zonoBundle_vertices
% test_zonoBundle_vertices - unit test function of vertices
%
% Syntax:
%    res = test_zonoBundle_vertices
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
V = vertices(zB);
res = isnumeric(V) && isempty(V) && size(V,1) == n;

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
V = vertices(zB);
V_true = [4 0; 1 3; -1 3; -2 2; -2 -1; 3 -1]';
res(end+1,1) = compareMatrices(V,V_true,1e-14);

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
V = vertices(zB);
res(end+1,1) = isnumeric(V) && isempty(V);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
