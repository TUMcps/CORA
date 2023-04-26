function res = test_zonoBundle_mptPolytope
% test_zonoBundle_mptPolytope - unit test function of mptPolytope
%
% Syntax:  
%    res = test_zonoBundle_mptPolytope
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
% Written:      24-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% fully-empty zonoBundle
zB = zonoBundle();
P = mptPolytope(zB);
res = isempty(P);

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
% convert to polytope
P = mptPolytope(zB);
% vertices
V = vertices(zB);
V_ = vertices(P);
% compare results
res(end+1,1) = compareMatrices(V,V_,1e-12);

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
% convert to polytope
P = mptPolytope(zB);
% vertices
V = vertices(zB);
V_ = vertices(P);
% compare results
res(end+1,1) = compareMatrices(V,V_,1e-12);

% combine results
res = all(res);

%------------- END OF CODE --------------
