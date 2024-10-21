function res = test_interval_length
% test_interval_length - unit_test_function of length
%
% Syntax:
%    res = test_interval_length
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
% Written:       19-January-2016
% Last update:   04-December-2023 (MW, add unbounded and matrix cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty
I = interval.empty(2);
assert(length(I) == 0);

% bounded
I = interval([-5, -4, -3, 0, 0, 5], [-2, 0, 2, 0, 5, 8]);
assert(length(I) == 6);
I = interval([-5; -4; -3; 0; 0; 5], [-2; 0; 2; 0; 5; 8]);
assert(length(I) == 6);

% bounded, matrix
I = interval([1 2 3; -2 1 -1],[3 4 6; -1 2 0]);
assert(length(I) == 3);

% unbounded
I = interval([-Inf;2],[1;Inf]);
assert(length(I) == 2);

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(length(I) == 3)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
