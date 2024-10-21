function res = test_polytope_zonoBundle
% test_polytope_zonoBundle - unit test function of zonoBundle
%
% Syntax:
%    res = test_polytope_zonoBundle
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       04-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-14;

% 2D, polytope is also a zonotope
A = [1 1; -1 1; -1 -1; 1 -1]; b = ones(4,1);
P = polytope(A,b);
V = vertices(P);
% convert to zonoBundle
zB = zonoBundle(P);
V_ = vertices(zB);
% compare vertices
assert(compareMatrices(V,V_,tol));


% 2D, bounded
A = [1 0; -1 1; -1 -1]; b = [0; 2; 0];
P = polytope(A,b);
V = vertices(P);
% convert to zonoBundle
zB = zonoBundle(P);
V_ = vertices(zB);
% compare vertices
assert(compareMatrices(V,V_,tol));


% 2D, fully empty -> unbounded
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
assertThrowsAs(@zonoBundle,'CORA:specialError',P);


% 2D, trivially fulfilled constraints -> unbounded
A = [0 0]; b = 1; Ae = [0 0]; be = 0;
P = polytope(A,b,Ae,be);
assertThrowsAs(@zonoBundle,'CORA:specialError',P);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
