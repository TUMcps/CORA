function res = test_interval_sqrt
% test_interval_sqrt - unit test function of sqrt
%
% Syntax:
%    res = test_interval_sqrt
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
% Written:       29-August-2019
% Last update:   04-December-2023 (MW, add empty and unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% empty
I = interval.empty(2);
I_sqrt = sqrt(I);
assert(representsa(I_sqrt,'emptySet'));

% bounded (only perfect squares)
I = interval([4; 9; 4; 16; 1], [9; 25; 36; 100; 4]);
I_sqrt = sqrt(I);
I_true = interval([2; 3; 2; 4; 1], [3; 5; 6; 10; 2]);
assert(isequal(I_sqrt,I_true,tol));

% unbounded
I = interval([2;4],[Inf;9]);
I_sqrt = sqrt(I);
I_true = interval([sqrt(2);2],[Inf;3]);
assert(isequal(I_sqrt,I_true,tol));

% out of bounds
I = interval(-2,1);
assertThrowsAs(@sqrt,'CORA:outOfDomain',I);

% n-d arrays
inf = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
sup = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = abs(interval(inf,sup));
Isqrt = sqrt(I);
inf = reshape([ 1.000 1.732 1.414 2.236 1.000 0.000 1.414 1.000 0.000 0.000 0.000 1.732 0.000 0.000 0.000 0.000 1.000 0.707 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
sup = reshape([ 1.225 2.000 2.000 3.162 1.732 0.000 1.732 1.414 1.000 1.414 1.414 2.000 0.000 0.000 0.000 0.000 1.414 1.000 1.732 1.414 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I_true = interval(inf,sup);
assert(isequal(Isqrt,I_true,1e-3))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
