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
res = true(0);

% empty
n = 2;
I = interval.empty(n);
lb = infimum(I);
res(end+1,1) = isempty(lb) && isnumeric(lb) && size(lb,1) == n;

% bounded
I = interval([-5, -4, -3, 0, 0, 5], [-2, 0, 2, 0, 5, 8]);
lb = infimum(I);
res(end+1,1) = all(withinTol(lb,[-5, -4, -3, 0, 0, 5],tol));

% unbounded
I = interval([-Inf;-2],[2;Inf]);
lb = infimum(I);
res(end+1,1) = all(withinTol(lb,[-Inf;-2],tol));

% matrix
I = interval([1 2; 3 4],[5 6; 7 8]);
lb = infimum(I);
res(end+1,1) = all(all(withinTol(lb,[1 2; 3 4],tol)));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
