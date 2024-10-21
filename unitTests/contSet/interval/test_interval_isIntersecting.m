function res = test_interval_isIntersecting
% test_interval_isIntersecting - unit test function of isIntersecting
%    note: only interval-to-interval tested!
%
% Syntax:
%    res = test_interval_isIntersecting
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
% Written:       12-March-2021
% Last update:   04-December-2023 (MW, add more cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
I1 = interval.empty(1);
I2 = interval(-1,1);
assert(~isIntersecting(I1,I2));
assert(~isIntersecting(I2,I1));

% bounded
I1 = interval([-2;-1],[1;2]);
I2 = interval([-4;-2],[-1;0]);
assert(isIntersecting(I1,I2));
assert(isIntersecting(I2,I1));
I1 = interval([-2;-1],[1;2]);
I2 = interval([-4;-2],[-3;0]);
assert(~isIntersecting(I1,I2));
assert(~isIntersecting(I2,I1));

% bounded and unbounded
I1 = interval([-2;-1],[1;2]);
I2 = interval([-4;1],[-1;Inf]);
assert(isIntersecting(I1,I2));
assert(isIntersecting(I2,I1));

% unbounded and unbounded
I1 = interval(-Inf,0);
I2 = interval(-1,Inf);
assert(isIntersecting(I1,I2));
assert(isIntersecting(I2,I1));

% check numeric
I = interval([2;4],[3;5]);
assert(all(isIntersecting(I,I.randPoint(10))));
assert(all(~isIntersecting(I,[1 1; 1 2])));

% check if tolerance is correctly used
I1 = interval([0;1],[1;2]);
I2 = interval([1;2],[2;3]) + 1e-10;
assert(isIntersecting(I1,I2));


% dimension mismatch
I1 = interval(-1,1);
I2 = interval([-1;-2],[2;1]);
assertThrowsAs(@isIntersecting,'CORA:dimensionMismatch',I1,I2);

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(isIntersecting(I,I))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
