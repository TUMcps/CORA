function res = test_zonotope_plus
% test_zonotope_plus - unit test function of plus
%
% Syntax:
%    res = test_zonotope_plus
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   09-August-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D zonotopes
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);
Z2 = zonotope([1 10; -1 -10]);

% Minkowski sum
Z_ = Z1 + Z2;

% obtain center and generator matrix
c_ = Z_.c; G_ = Z_.G;

% compare to true result
true_c = [-3; 0];
true_G = [-3, -2, -1, 10; ...
            2, 3, 4, -10];
res(end+1,1) = compareMatrices(c_,true_c) && compareMatrices(G_,true_G);

% Minkowski sum with empty set
Z_empty = zonotope.empty(2);
res(end+1,1) = representsa(Z1 + Z_empty,'emptySet');


% add results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
