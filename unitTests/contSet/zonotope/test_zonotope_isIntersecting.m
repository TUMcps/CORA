function res = test_zonotope_isIntersecting
% test_zonotope_isIntersecting - unit test function of isIntersecting
%
% Syntax:
%    res = test_zonotope_isIntersecting
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

% Authors:       Tobias Ladner
% Written:       20-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create zonotope
c = [-4; 1];
G = [-3, -2, -1; 2, 3, 4];
Z1 = zonotope(c,G);

c = [-8;4];
G = [1 0 1; 0 1 1];
Z2 = zonotope(c,G);

% should intersect
assert(isIntersecting(Z1,Z2));

% move outside
Z2 = Z2 - [4;5];
assert(~isIntersecting(Z1,Z2));

% check containment 'is intersecting'
Z2 = enlarge(Z1,[2;1.5]);
assert(isIntersecting(Z1,Z2));

% check touching zonotopes
c = [8; 1];
G = [3, 2, 1; 2, 3, 4];
Z2 = zonotope(c,G) + 1e-10;
assert(isIntersecting(Z1,Z2));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
