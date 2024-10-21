function res = test_interval_center
% test_interval_center - unit test function of center
%
% Syntax:
%    res = test_interval_center
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

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       15-January-2016
% Last update:   03-December-2023 (MW, add unbounded case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% 1. empty set
n = 2;
I = interval.empty(n);
c = center(I);
assert(isempty(c) && isnumeric(c) && all(size(c) == [2,0]));

% 2. bounded
I = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = center(I);
c_true = [-3.5,-2,-0.5,0,2.5,6.5];
assert(all(withinTol(c,c_true,tol)));

% 3. unbounded
I = interval(-Inf,2);
c = center(I);
assert(isnan(c));
I = interval(2,Inf);
c = center(I);
assert(isnan(c));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
c = center(I);
c_true = (lb+ub)/2;
assert(all(withinTol(c,c_true),"all"))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
