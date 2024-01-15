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
res = true(0);

% empty interval
n = 2;
I = interval.empty(n);
r = rad(I);
res(end+1,1) = isempty(r) && isnumeric(r) && size(r,1) == n;

% 6D interval
I = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
r = rad(I);
r_true = [1.5, 2, 2.5, 0, 2.5, 1.5];
res(end+1,1) = all(withinTol(r,r_true,tol));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
