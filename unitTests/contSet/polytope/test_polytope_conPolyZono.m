function res = test_polytope_conPolyZono
% test_polytope_conPolyZono - unit test function of conversion to
%    constrained polynomial zonotopes
%
% Syntax:
%    res = test_polytope_conPolyZono
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

% init polytope
A = [2 1; -1 3; -2 -2; 1 -3];
b = ones(4,1);
P = polytope(A,b);

% convert to polyZonotope
cPZ = conPolyZono(P);

% sample random points from converted conPolyZono and check if all are
% contained in the polytope
p_pZ = randPoint(cPZ,100);

res = all(contains(P,p_pZ));

% ------------------------------ END OF CODE ------------------------------
