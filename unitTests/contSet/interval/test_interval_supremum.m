function res = test_interval_supremum
% test_interval_supremum - unit test function of the supremum of an interval
%
% Syntax:
%    res = test_interval_supremum
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
% Last update:   04-December-2023 (MW, add empty and unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define problem
tol = 1e-9;

% empty
n = 2;
I = interval.empty(n);
ub = supremum(I);
assert(isnumeric(ub) && isempty(ub) && size(ub,1) == n);

% bounded
I = interval([-5, -4, -3, 0, 0, 5], [-2, 0, 2, 0, 5, 8]);
ub = supremum(I);
assert(all(withinTol(ub,[-2, 0, 2, 0, 5, 8],tol)));

% unbounded
I = interval([-Inf;-2],[2;Inf]);
ub = supremum(I);
assert(all(withinTol(ub,[2;Inf],tol)));

% matrix
I = interval([1 2; 3 4],[5 6; 7 8]);
ub = supremum(I);
assert(all(all(withinTol(ub,[5 6; 7 8],tol))));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(all(withinTol(I.sup,ub),'all'))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
