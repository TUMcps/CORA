function res = test_zonotope_quadMap
% test_zonotope_quadMap - unit test function of quadMap
%
% Syntax:
%    res = test_zonotope_quadMap
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

% Authors:       Matthias Althoff
% Written:       26-July-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create zonotopes
Z1 = zonotope([-4, -3, -2; 1, 2, 3]);
Z2 = zonotope([1, 4, 2; -3, 2, -1]);

% create matrices
Q{1} = [1 -2; -3 4];
Q{2} = [0.5 0; 2 -1];

% 1. quadMapSingle: Z1*Q*Z1 -----------------------------------------------

% obtain result
Zres = quadMap(Z1,Q);

% obtain center and generator matrix
c = Zres.c;
G = Zres.G;

% true result
true_c = [102.5; -16.25];
true_G = [27.5, 35, 95, 110, 125; ...
            -5.75, -9.5, -14, -26, -32];

% compare solutions
res(1) = compareMatrices(c,true_c) && compareMatrices(G,true_G);

% 2. quadMapMixed: Z1*Q*Z2 ------------------------------------------------

% obtain result
Zres = quadMap(Z1,Z2,Q);

% obtain center and generator matrix
c = Zres.c;
G = Zres.G;

% true result
true_c = [-43; 3];
true_G = [-51, -59, -4, -8, -12, -26, -32, -38;
             8.5, 14, -2, 6, 14, 1, 7, 13];

% compare solutions
res(2) = compareMatrices(c,true_c) && compareMatrices(G,true_G);


% gather results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
