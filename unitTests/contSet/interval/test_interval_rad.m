function res = test_interval_rad
% test_interval_rad - unit_test_function of rad
%
% Syntax:
%    res = test_interval_rad
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

% Authors:       Dmitry Grebenyuk
% Written:       19-January-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% empty interval
n = 2;
I = interval.empty(n);
r = rad(I);
assert(isempty(r) && isnumeric(r) && size(r,1) == n);

% 6D interval
I = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
r = rad(I);
r_true = [1.5, 2, 2.5, 0, 2.5, 1.5];
assert(all(withinTol(r,r_true,tol)));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
r_true = reshape([ 0.250 0.500 1.000 2.500 1.000 0.000 0.500 0.500 0.500 1.000 1.500 0.500 0.000 0.000 0.000 0.000 0.500 0.250 1.000 1.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
assert(isequal(rad(I),r_true))


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
