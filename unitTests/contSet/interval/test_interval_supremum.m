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
res = true(0);

% empty
n = 2;
I = interval.empty(n);
ub = supremum(I);
res(end+1,1) = isnumeric(ub) && isempty(ub) && size(ub,1) == n;

% bounded
I = interval([-5, -4, -3, 0, 0, 5], [-2, 0, 2, 0, 5, 8]);
ub = supremum(I);
res(end+1,1) = all(withinTol(ub,[-2, 0, 2, 0, 5, 8],tol));

% unbounded
I = interval([-Inf;-2],[2;Inf]);
ub = supremum(I);
res(end+1,1) = all(withinTol(ub,[2;Inf],tol));

% matrix
I = interval([1 2; 3 4],[5 6; 7 8]);
ub = supremum(I);
res(end+1,1) = all(all(withinTol(ub,[5 6; 7 8],tol)));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
