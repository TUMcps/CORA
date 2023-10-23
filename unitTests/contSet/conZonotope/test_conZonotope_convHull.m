function res = test_conZonotope_convHull
% test_conZonotope_convHull - unit test function for convex hull of a
%    constrained zonotope with other sets
%
% Syntax:
%    res = test_conZonotope_convHull
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
% Written:       21-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% instantiate constrained zonotopes
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1]; b = 1;
cZ1 = conZonotope(Z,A,b);

Z = [4 2 0 0;4 1 1 0];
A = [1 1 -1]; b = 0;
cZ2 = conZonotope(Z,A,b);

% compute convex hull
cZ_ = convHull(cZ1,cZ2);

% compute vertices
V1 = vertices(cZ1);
V2 = vertices(cZ2);
V_ = vertices(cZ_);

% check vertices
resvec(end+1) = compareMatrices(V_,[V1,V2],1e-14,'subset');

% convert second constrained zonotope to polytope and check again
P2 = polytope(cZ2);
cZ_ = convHull(cZ1,P2);


% check vertices
V_ = vertices(cZ_);
resvec(end+1) = compareMatrices(V_,[V1,V2],1e-12,'subset');

% add a third constrained zonotope
Z = [-4 3 0 1; 5 0 2 1];
A = [1 0 1]; b = 1;
cZ3 = conZonotope(Z,A,b);

% compute convex hull
cZ_ = convHull(cZ1,{cZ2,cZ3});

% compute vertices
V3 = vertices(cZ3);
V_ = vertices(cZ_);

% check vertices
resvec(end+1) = compareMatrices(V_,[V1,V2,V3],1e-14,'subset');

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
