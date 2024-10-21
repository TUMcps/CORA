function res = test_ellipsoid_norm
% test_ellipsoid_norm - unit test function of norm
%
% Syntax:
%    res = test_ellipsoid_norm
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
 
% assume true
res = true;

% init cases
E1 = ellipsoid([ 5.4387811500952807 12.4977183618314545 ; 12.4977183618314545 29.6662117284481646 ], [ -0.7445068341257537 ; 3.5800647524843665 ], 0.000001);
Ed1 = ellipsoid([ 4.2533342807136076 0.6346400221575308 ; 0.6346400221575309 0.0946946398147988 ], [ -2.4653656883489115 ; 0.2717868749873985 ], 0.000001);
E0 = ellipsoid([ 0.0000000000000000 0.0000000000000000 ; 0.0000000000000000 0.0000000000000000 ], [ 1.0986933635979599 ; -1.9884387759871638 ], 0.000001);

% remove center
E1 = E1 - center(E1);
Ed1 = Ed1 - center(Ed1);
E0 = E0 - center(E0);

% empty set
assert(norm(ellipsoid.empty(2)) == -Inf);

TOL = E1.TOL;
n = dim(E1);
N = 10*n;
norm1 = max(sqrt(sum(randPoint(E1,N,'extreme').^2,1)));
normd = max(sqrt(sum(randPoint(Ed1,N,'extreme').^2,1)));
norm0 = max(sqrt(sum(randPoint(E0,N,'extreme').^2,1)));

assert(norm(E1)+TOL >= norm1)
assert(norm(Ed1)+TOL >= normd)
assert(norm(E0)+TOL >= norm0)

end

% ------------------------------ END OF CODE ------------------------------
