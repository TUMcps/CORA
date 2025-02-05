function res = test_zonotope_contains
% test_zonotope_contains - unit test function of contains
%
% Syntax:
%    res = test_zonotope_contains
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
% Written:       12-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% point-in-zono ---

% 2D zonotopes
c = [0; 1]; G = [1 2 1; -1 0 1];
Z = zonotope(c,G);
assert(contains(Z, [1.5797;1.3568]))
assert(~contains(Z, [3.5797;2.3568]))

% degenerate
Z = zonotope([-1;0], [1;1]);
assert(contains(Z, [-0.5;0.5]))
assert(~contains(Z, [-0.5;-0.5]))

% almost degenerate
Z = zonotope([-1;0], [1 1e-5;1 2e-5]);
assert(contains(Z, [-0.5;0.5]))
assert(~contains(Z, [-0.5;-0.5]))

% zono-in-zono ---

% empty zonotope
% Z = zonotope.empty(2);
% res(end+1,1) = contains(Z,Z);

% 2D zonotopes
c = [0; 1]; G = [1 2 1; -1 0 1];
Z1 = zonotope(c,G);
c = [-1; 1.5]; G = [0.2 0; -0.1 0.1];
Z2 = zonotope(c,G);
assert(contains(Z1,Z2));
assert(~contains(Z2,Z1));

% inner zonotope is just a point
Z2 = zonotope(c);
assert(contains(Z1,Z2));

% choose LP method for containment
assert(contains(Z1,Z2,'approx:st'));

% both zonotope are just points
Z1 = zonotope(zeros(4,1));
Z2 = zonotope(ones(4,1));
assert(contains(Z1,Z1));
assert(~contains(Z1,Z2));

% two sets with one being degenerate
c = [ 5.000 ; 0.000 ];
G = [ 0.000 ; 1.000 ];
Z1 = zonotope(c,G);
c = [ 5.650 ; 0.000 ];
G = [ 0.000 0.050 0.000 0.000 0.000 ; 0.937 0.000 -0.005 -0.000 0.000 ];
Z2 = zonotope(c,G);
assert(~contains(Z1,Z2))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
