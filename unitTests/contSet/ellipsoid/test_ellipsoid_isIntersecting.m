function res = test_ellipsoid_isIntersecting
% test_ellipsoid_isIntersecting - unit test function of isIntersecting
%
% Syntax:
%    res = test_ellipsoid_isIntersecting
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
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init cases
E1 = ellipsoid([ 5.4387811500952807 12.4977183618314545 ; 12.4977183618314545 29.6662117284481646 ], [ -0.7445068341257537 ; 3.5800647524843665 ], 0.000001);
E0 = ellipsoid([ 0.0000000000000000 0.0000000000000000 ; 0.0000000000000000 0.0000000000000000 ], [ 1.0986933635979599 ; -1.9884387759871638 ], 0.000001);

res = true;

% empty set: rewrite using emptySet class
assert(~isIntersecting(E1,ellipsoid.empty(2)));

E0 = ellipsoid(E0.Q,E1.q);
assert(isIntersecting(E1,E0))

r1 = rad(interval(E1));
E2 = ellipsoid(E1.Q,E1.q+2*r1);
assert(~isIntersecting(E2,E1))
    
end

% ------------------------------ END OF CODE ------------------------------
