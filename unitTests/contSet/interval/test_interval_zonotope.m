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

res = true(0);
tol = 1e-9;

% empty interval
I = interval.empty(2);
Z = zonotope(I);
res(end+1,1) = representsa(Z,'emptySet') && dim(Z) == 2;

% bounded interval
I = interval([-2;-5],[3;1]);
Z = zonotope(I);
c_true = [0.5;-2]; G_true = [2.5 0; 0 3];
res(end+1,1) = all(withinTol(center(Z),c_true,tol));
res(end+1,1) = compareMatrices(generators(Z),G_true,tol);

% unbounded interval
I = interval([-Inf;-2],[4;2]);
try
    zonotope(I);
    res = false;
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
