function res = test_ellipsoid_randPoint
% test_ellipsoid_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_ellipsoid_randPoint
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
E2 = ellipsoid([ 0.3542431574242590 0.0233699257103926 ; 0.0233699257103926 2.4999614009532856 ], [ 0.0873801375346114 ; -2.4641617305288825 ], 0.000001);
E0 = ellipsoid([ 0.0000000000000000 0.0000000000000000 ; 0.0000000000000000 0.0000000000000000 ], [ 1.0986933635979599 ; -1.9884387759871638 ], 0.000001);

% empty set
n = 2;
E = ellipsoid.empty(n);
p = randPoint(E);
assert(isempty(p) && isnumeric(p) && size(p,1) == n);


n = dim(E1);
N = 5*n;

assert(aux_withinRange(E1,N))
assert(aux_withinRange(Ed1,N))
assert(aux_withinRange(E0,N))

end


% Auxiliary functions -----------------------------------------------------

function res = aux_withinRange(E,N)
    % all extreme points need to be between min(radius) and max(radius)
    n = dim(E);
    Y = randPoint(E,N,'extreme');
    nY = sqrt(sum((Y-E.q).^2,1));
    rE = radius(E,n);
    IntR = interval(min(rE),max(rE));
    res = all(contains(IntR,nY));
end

% ------------------------------ END OF CODE ------------------------------
