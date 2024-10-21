function res = test_ellipsoid_or
% test_ellipsoid_or - unit test function of or
%
% Syntax:
%    res = test_ellipsoid_or
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
E2 = ellipsoid([ 0.3542431574242590 0.0233699257103926 ; 0.0233699257103926 2.4999614009532856 ], [ 0.0873801375346114 ; -2.4641617305288825 ], 0.000001);
E0 = ellipsoid([ 0.0000000000000000 0.0000000000000000 ; 0.0000000000000000 0.0000000000000000 ], [ 1.0986933635979599 ; -1.9884387759871638 ], 0.000001);


% empty set
E_empty = ellipsoid.empty(2);
assert(isequal(or(E1,E_empty),E1))

% test non-deg
Eres_nd = E1 | E2;
Y_nd = [randPoint(E1,2),randPoint(E2,2)];
assert(all(contains(Eres_nd,Y_nd)))

% test zero rank ellipsoid
Eres_0 = E1 | E0;
Y_0 = [randPoint(E1,2),E0.q];
assert(all(contains(Eres_0,Y_0)))
    
end

% ------------------------------ END OF CODE ------------------------------
