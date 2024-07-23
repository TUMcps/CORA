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

res = true(0);
tol = 1e-12;

% empty case
P = polytope.empty(2);
Z = zonotope(P);
res(end+1,1) = representsa(Z,'emptySet') && dim(Z) == 2;

% 1D, bounded
A = [1; -1]; b = [2; -1];
P = polytope(A,b);
Z = zonotope(P);
Z_true = zonotope(1.5,0.5);
res(end+1,1) = Z == Z_true;

% 2D, bounded
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
Z = zonotope(P);
res(end+1,1) = contains(Z,P,'exact',tol);

% 2D, bounded, vertex instantiation
V = [-1 0; 1 2; 1 -2]';
P = polytope(V);
Z = zonotope(P);
res(end+1,1) = contains(Z,P,'exact',tol);


% combine results
res = all(res);


% unsupported
try
    % unbounded
    P = polytope.Inf(2);
    Z = zonotope(P);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
