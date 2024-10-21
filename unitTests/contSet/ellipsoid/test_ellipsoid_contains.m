function res = test_ellipsoid_contains
% test_ellipsoid_contains - unit test function of containment checks
%
% Syntax:
%    res = test_ellipsoid_contains
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

% Authors:       Victor Gassmann, Mark Wetzlinger
% Written:       26-July-2021
% Last update:   ---
% Last revision: 22-September-2024 (MW, explicit cases)

% ------------------------------ BEGIN CODE -------------------------------

%%% all-zero shape matrix
Q = zeros(3,3); q = [3;2;1];
E = ellipsoid(Q,q);

% empty set
O = emptySet(3);
assert(contains(E,O));

% fullspace
fs = fullspace(3);
assert(~contains(E,fs));

% points
assert(contains(E,q));
assert(~contains(E,q+q));

% zonotope
Z = zonotope(q);
assert(contains(E,Z));

% polytope
P = polytope(q);
assert(contains(E,P));


%%% non-degenerate ellipsoid
Q = [2 0; 0 1]; q = [2;-1];
M = [2 -1; 1 1];
E = M * ellipsoid(Q,q);

% points
p = [4 0; 3 1; 5 1; 6 2; 8 2]';
assert(all(contains(E,p)));

% zonotope
c = [5;1]; G_inside = [0.2 0.3 0.5; 0 -0.3 0.2]; G_outside = [1 1 0.5; 0 -1 1];
Z_inside = zonotope(c,G_inside);
assert(contains(E,Z_inside));
Z_outside = zonotope(c,G_outside); 
assert(~contains(E,Z_outside));

% interval
I_inside = interval([3;0.5],[6;1]);
assert(contains(E,I_inside));
I_outside = interval([3;0],[7;2]);
assert(~contains(E,I_outside));

% polytope
P_inside = polytope([1 0; -1 1; -1 -1],[4;-2.2;-4]);
assert(contains(E,P_inside));
P_outside = polytope([1 0; -1 1; -1 -1],[7;-2;-4]);
assert(~contains(E,P_outside));


%%% degenerate ellipsoid
Q = [2 0 0; 0 3 0; 0 0 1]; q = [2;-1;4];
M = [2 -1 1; 1 1 2; -1 3 2];
E = M * ellipsoid(Q,q);

% points
p = M * (q + [0.5 0.2 -0.5; -0.8 -1.2 0.2; 1.2 1.5 -0.5]');
assert(all(contains(E,p)));

% zonotope
Z = M * zonotope(q,[0.1 0.2 0.3; -0.2 0.2 0.1; 0.0 0.1 -0.1]);
assert(contains(E,Z));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
