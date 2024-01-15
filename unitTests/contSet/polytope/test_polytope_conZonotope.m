function res = test_polytope_conZonotope
% test_polytope_conZonotope - unit test function of conversion to
%    constrained zonotopes
%
% Syntax:
%    res = test_polytope_conZonotope
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
% Written:       29-November-2022
% Last update:   14-July-2023 (MW, add empty cases)
%                27-July-2023 (MW, add more cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% fully empty polytope
P = polytope.empty(2);
cZ = conZonotope(P);
res(end+1,1) = isemptyobject(cZ) && dim(cZ) == 2;

% 1D, only inequalities, bounded
A = [1;-1]; b = [2;5];
P = polytope(A,b);
V = vertices(P);
cZ = conZonotope(P);
V_ = vertices(cZ);
res(end+1,1) = compareMatrices(V,V_,1e-14);

% 1D, single point
Ae = 1; be = 4;
P = polytope([],[],Ae,be);
V = vertices(P);
cZ = conZonotope(P);
V_ = vertices(cZ);
res(end+1,1) = compareMatrices(V,V_,1e-14);


% 2D, only inequalities, empty
A = [1 1; -1 1; -1 -1; 1 -1]; b = [1; 1; 1; -2];
P = polytope(A,b);
% convert to constrained zonotope
cZ = conZonotope(P);
res(end+1,1) = representsa(cZ,'emptySet');
cZ = conZonotope(P,'exact:vertices');
res(end+1,1) = representsa(cZ,'emptySet');

% 2D, bounded, vertex instantiation
V = [3 -2; 3 2; 1 3; 1 -1]';
P = polytope(V);
% convert to constraiend zonotope with default method
cZ = conZonotope(P);
% alternative solution
% Z = [0 3 0 1;0 0 2 1]; A = [1 0 1]; b = 1;
% cZ = conZonotope(Z,A,b);
V_ = vertices(cZ);
res(end+1,1) = compareMatrices(V,V_,1e-14);
% convert to constrained zonotope with special method
cZ_vert = conZonotope(P,'exact:vertices');
V_vert = vertices(cZ_vert);
res(end+1,1) = compareMatrices(V,V_vert,1e-14);

% 2D, only inequalities
A = [2 1; -1 1; -1 -3; 4 -1];
A = (A' ./ vecnorm(A'))'; b = ones(4,1);
P = polytope(A,b);
% compute vertices
V = vertices(P);
% convert to constrained zonotope with default method
cZ = conZonotope(P);
V_ = vertices(cZ);
res(end+1,1) = compareMatrices(V,V_,1e-8);
% convert to constrained zonotope with special method
cZ_vert = conZonotope(P,'exact:vertices');
V_vert = vertices(cZ_vert);
res(end+1,1) = compareMatrices(V,V_vert,1e-8);

% 2D, degenerate, inequalities and equalities
A = [0 1; 0 -1]; b = [3;-1]; Ae = [1 1]; be = 1;
P = polytope(A,b,Ae,be);
V = vertices(P);
% convert to constrained zonotope with default method
cZ = conZonotope(P);
V_ = vertices(cZ);
res(end+1,1) = compareMatrices(V,V_,1e-14);
% convert to constrained zonotope with special method
cZ_vert = conZonotope(P,'exact:vertices');
V_vert = vertices(cZ_vert);
res(end+1,1) = compareMatrices(V,V_vert,1e-14);

% combine results
res = all(res);


% impossible conversions

% 2D, fully empty (fullspace)
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
try
    cZ = conZonotope(P);
    res = false;
end

% 2D, unbounded
A = [1 0]; b = 1;
P = polytope(A,b);
try
    cZ = conZonotope(P);
    res = false;
end
try
    cZ = conZonotope(P,'exact:vertices');
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
