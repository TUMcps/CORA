function res = test_interval_isscalar
% test_interval_isscalar - unit test function of isscalar
%
% Syntax:
%    res = test_interval_isscalar
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
% Written:       16-January-2016
% Last update:   04-December-2023 (MW, add unbounded and matrix cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty
I = interval.empty(1);
assert(~isscalar(I));

% bounded, scalar
I = interval(-5, 2);
assert(isscalar(I));

% bounded, vector
I = interval([-5,-4,-3,0,0,5],[-2,0,2,0,5,8]);
assert(~isscalar(I));

% bounded, matrix
I = interval([1 2; 3 4],[5 6; 7 8]);
assert(~isscalar(I));

% unbounded, scalar
I = interval(-Inf,Inf);
assert(isscalar(I));

% unbounded, vector
I = interval([-Inf;1],[2;Inf]);
assert(~isscalar(I));

% unbounded, matrix
I = interval([-Inf 2; 1 -Inf],[2 4; Inf 0]);
assert(~isscalar(I));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(~isscalar(I));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
