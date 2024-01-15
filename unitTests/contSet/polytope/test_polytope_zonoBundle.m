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

res = true(0);

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
res(end+1,1) = compareMatrices(V,V_,tol);


% 2D, bounded
A = [1 0; -1 1; -1 -1]; b = [0; 2; 0];
P = polytope(A,b);
V = vertices(P);
% convert to zonoBundle
zB = zonoBundle(P);
V_ = vertices(zB);
% compare vertices
res(end+1,1) = compareMatrices(V,V_,tol);


% combine results
res = all(res);


% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
try
    zB = zonoBundle(P);
    res = false;
end

% 2D, trivially fulfilled constraints
A = [0 0]; b = 1; Ae = [0 0]; be = 0;
P = polytope(A,b,Ae,be);
try
    zB = zonoBundle(P);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
