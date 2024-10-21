function res = test_interval_isequal
% test_interval_isequal - unit test function of isequal
%
% Syntax:
%    res = test_interval_isequal
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   03-December-2022 (MW, add Inf case)
%                23-December-2022 (MW, add matrix case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% bounded, different size
I1 = interval([3;2;1],[4;5;6]);
I2 = interval([3;2],[4;5]);
assert(~isequal(I1,I2));

% bounded, same size
I1 = interval([-3; -9; -4; -7; -1],[4; 2; 6; 3; 8]);
I2 = interval([-3; -9; -4; -7; -1],[4; 2; 6; 2; 8]);
assert(isequal(I1,I1));
assert(~isequal(I1,I2));

% unbounded
I1 = interval([-Inf; -9; -4; -7; -1], [4; 2; 6; Inf; 8]);
I2 = interval([-Inf; -9; -4; -7; -1], [4; 2; Inf; 2; 8]);
assert(isequal(I1,I1));
assert(~isequal(I1,I2));

% instantiate matrix interval
I1 = interval([-2 -3 -Inf; -4 -Inf -1],[2 Inf 5; Inf 3 0]);
I2 = interval([-2 -3 -Inf; -4 -Inf -1],[2 Inf Inf; Inf 3 0]);
I3 = I1([1,2],[1,2]);
assert(isequal(I1,I1));
assert(~isequal(I1,I2));
assert(~isequal(I1,I3));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(isequal(I,I))
assert(~isequal(I,I+1))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
