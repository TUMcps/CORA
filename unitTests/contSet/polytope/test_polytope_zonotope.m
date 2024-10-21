function res = test_polytope_zonotope
% test_polytope_zonotope - unit test function of zonotope conversion
%
% Syntax:
%    res = test_polytope_zonotope
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
% Written:       08-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-12;

% empty case
P = polytope.empty(2);
Z = zonotope(P);
assert(representsa(Z,'emptySet') && dim(Z) == 2);

% 1D, bounded
A = [1; -1]; b = [2; -1];
P = polytope(A,b);
Z = zonotope(P);
Z_true = zonotope(1.5,0.5);
assert(Z == Z_true);

% 2D, bounded
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
Z = zonotope(P);
assert(contains(Z,P,'exact',tol));

% 2D, bounded, vertex instantiation
V = [-1 0; 1 2; 1 -2]';
P = polytope(V);
Z = zonotope(P);
assert(contains(Z,P,'exact',tol));


% unsupported: unbounded
P = polytope.Inf(2);
assertThrowsAs(@zonotope,'CORA:specialError',P);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
