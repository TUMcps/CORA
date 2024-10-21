function res = test_ellipsoid_enclose
% test_ellipsoid_enclose - unit test function of enclose
%
% Syntax:
%    res = test_ellipsoid_enclose
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% init cases
E1 = ellipsoid([ 5.4387811500952807 12.4977183618314545 ; 12.4977183618314545 29.6662117284481646 ], [ -0.7445068341257537 ; 3.5800647524843665 ], 0.000001);
Ed1 = ellipsoid([ 4.2533342807136076 0.6346400221575308 ; 0.6346400221575309 0.0946946398147988 ], [ -2.4653656883489115 ; 0.2717868749873985 ], 0.000001);
E2 = ellipsoid([ 0.3542431574242590 0.0233699257103926 ; 0.0233699257103926 2.4999614009532856 ], [ 0.0873801375346114 ; -2.4641617305288825 ], 0.000001);
E0 = ellipsoid([ 0.0000000000000000 0.0000000000000000 ; 0.0000000000000000 0.0000000000000000 ], [ 1.0986933635979599 ; -1.9884387759871638 ], 0.000001);
n = length(E1.q);

% test non-deg ellipsoid
E = enclose(E1,E2);
Y = [randPoint(E1,2*n),randPoint(E2,2*n)];
assert(all(contains(E,Y)))

% test degenerate ellipsoids
Ed = enclose(Ed1,E0);
Yd = [randPoint(Ed1,2*n),E0.q];
assert(all(contains(Ed,Yd)))
    
end

% ------------------------------ END OF CODE ------------------------------
