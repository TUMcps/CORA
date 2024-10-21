function res = test_interval_infimum
% test_interval_infimum - unit test function of the infimum of an interval
%
% Syntax:
%    res = test_interval_infimum
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

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       14-January-2016
% Last update:   04-December-2023
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define problem
tol = 1e-9;

% empty
n = 2;
I = interval.empty(n);
lb = infimum(I);
assert(isempty(lb) && isnumeric(lb) && size(lb,1) == n);

% bounded
I = interval([-5, -4, -3, 0, 0, 5], [-2, 0, 2, 0, 5, 8]);
lb = infimum(I);
assert(all(withinTol(lb,[-5, -4, -3, 0, 0, 5],tol)));

% unbounded
I = interval([-Inf;-2],[2;Inf]);
lb = infimum(I);
assert(all(withinTol(lb,[-Inf;-2],tol)));

% matrix
I = interval([1 2; 3 4],[5 6; 7 8]);
lb = infimum(I);
assert(all(all(withinTol(lb,[1 2; 3 4],tol))));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(all(withinTol(I.inf,lb),'all'))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
