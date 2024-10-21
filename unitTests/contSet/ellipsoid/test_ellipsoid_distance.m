function res = test_ellipsoid_distance
% test_ellipsoid_distance - unit test function of distance
%
% Syntax:
%    res = test_ellipsoid_distance
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

% Authors:       Victor Gassmann
% Written:       26-July-2021
% Last update:   07-July-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
% init cases
E1 = ellipsoid([ 5.4387811500952807 12.4977183618314545 ; 12.4977183618314545 29.6662117284481646 ], [ -0.7445068341257537 ; 3.5800647524843665 ], 0.000001);
E0 = ellipsoid([ 0.0000000000000000 0.0000000000000000 ; 0.0000000000000000 0.0000000000000000 ], [ 1.0986933635979599 ; -1.9884387759871638 ], 0.000001);


% check all-zero ellipsoid
assert(withinTol(distance(E1,E0),distance(E1,E0.q),E1.TOL))

n = length(E1.q);
% check polytope: construct hyperplane through E1.q
l1 = randn(1,n);
H = polytope([],[],l1,l1*E1.q);
% check if distance is <0
assert(distance(E1,H) < 0)

% check polytope: construct second hyperplane (also contains center
% of ellipsoid)
l2 = randn(1,n);
P = polytope([l1;l2],[l1*E1.q;l2*E1.q]);
assert(distance(E1,P) <= 1e-6)

end

% ------------------------------ END OF CODE ------------------------------
