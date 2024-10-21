function res = test_ellipsoid_cartProd
% test_ellipsoid_cartProd - unit test function of cartProd
%
% Syntax:
%    res = test_ellipsoid_cartProd
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
% Last update:   03-January-2023 (MW, add ellipsoid-numeric case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
% init cases
E1 = ellipsoid([ 5.4387811500952807 12.4977183618314545 ; 12.4977183618314545 29.6662117284481646 ], [ -0.7445068341257537 ; 3.5800647524843665 ], 0.000001);
Ed1 = ellipsoid([ 4.2533342807136076 0.6346400221575308 ; 0.6346400221575309 0.0946946398147988 ], [ -2.4653656883489115 ; 0.2717868749873985 ], 0.000001);
E2 = ellipsoid([ 0.3542431574242590 0.0233699257103926 ; 0.0233699257103926 2.4999614009532856 ], [ 0.0873801375346114 ; -2.4641617305288825 ], 0.000001);
E0 = ellipsoid([ 0.0000000000000000 0.0000000000000000 ; 0.0000000000000000 0.0000000000000000 ], [ 1.0986933635979599 ; -1.9884387759871638 ], 0.000001);

% empty set
% try 
%     cartProd(E_c{1}.E1,ellipsoid.empty(2));
%     res = false;
% catch ME 
%     if ~strcmp(ME.identifier,'CORA:notSupported')
%         rethrow(ME);
%     else
%         res = true;
%     end
% end

% test non-deg
Eres_nd = cartProd(E1,E2);
Y_nd = combineVec(randPoint(E1,2),randPoint(E2,2));
assert(all(contains(Eres_nd,Y_nd)))

% test deg
Eres_d = cartProd(E1,Ed1);
Y_d = combineVec(randPoint(E1,2),randPoint(Ed1,2));
assert(all(contains(Eres_d,Y_d)))

% test zero rank ellipsoid
Eres_0 = cartProd(E1,E0);
Y_0 = combineVec(randPoint(E1,2),randPoint(E0,2));
assert(all(contains(Eres_0,Y_0)))

% ellipsoid-numeric case
E = ellipsoid(eye(2),[1;2]);
num = -1;

% compute Cartesian product
E_ = cartProd(E,num);
E__ = cartProd(num,E);

% check result
assert(compareMatrices(E_.Q,[1 0 0; 0 1 0; 0 0 0]))
assert(compareMatrices(E_.q,[1;2;-1]))
assert(compareMatrices(E__.Q,[0 0 0; 0 1 0; 0 0 1]))
assert(compareMatrices(E__.q,[-1;1;2]))

end

% ------------------------------ END OF CODE ------------------------------
