function res = test_polytope_polyZonotope
% test_polytope_polyZonotope - unit test function of conversion to
%    polynomial zonotopes
%
% Syntax:
%    res = test_polytope_polyZonotope
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

% 2D, bounded
A = [2 1; -1 3; -2 -2; 1 -3]; b = ones(4,1);
P = polytope(A,b);
% convert to polyZonotope
pZ = polyZonotope(P);
% sample random points from converted polyZonotope and check if all are
% contained in the polytope
p_pZ = randPoint(pZ,1000);
res(end+1,1) = all(contains(P,p_pZ));


% combine results
res = all(res);


% 3D, fully empty
A = zeros(0,3); b = zeros(0,0);
P = polytope(A,b);
try
    pZ = polyZonotope(P);
    res = false;
end

% 2D, trivially fulfilled constraints
A = [0 0]; b = 1; Ae = [0 0]; be = 0;
P = polytope(A,b,Ae,be);
try
    pZ = polyZonotope(P);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
