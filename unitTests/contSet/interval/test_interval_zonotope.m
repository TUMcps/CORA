function res = test_interval_zonotope
% test_interval_zonotope - unit test function of zonotope conversion
%
% Syntax:
%    res = test_interval_zonotope
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
% Written:       28-August-2019
% Last update:   04-December-2023 (MW, add empty and unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-9;

% empty interval
I = interval.empty(2);
Z = zonotope(I);
assert(representsa(Z,'emptySet') && dim(Z) == 2);

% bounded interval
I = interval([-2;-5],[3;1]);
Z = zonotope(I);
c_true = [0.5;-2]; G_true = [2.5 0; 0 3];
assert(all(withinTol(center(Z),c_true,tol)));
assert(compareMatrices(generators(Z),G_true,tol));

% unbounded interval
I = interval([-Inf;-2],[4;2]);
assertThrowsAs(@zonotope,'CORA:wrongValue',I);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
