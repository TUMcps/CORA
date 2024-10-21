function res = test_interval_eq
% test_interval_eq - unit test function of eq, overloaded '==' operator
%
% Syntax:
%    res = test_interval_eq
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
% Written:       29-August-2019
% Last update:   04-December-2023 (MW, add empty and unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
I1 = interval.empty(2);
assert(I1 == I1);


% bounded
I1 = interval([-2; -4; -3],[2; 3; 1]);
I2 = interval([-3; 0; -4],[2; 3; 1]);
assert(I1 == I1);
assert(~(I1 == I2));

% unbounded
I1 = interval(-Inf,0);
I2 = interval(-Inf,Inf);
I3 = interval(0,Inf);
assert(I1 == I1);
assert(I2 == I2);
assert(I3 == I3);
assert(~(I1 == I2));
assert(~(I1 == I3));
assert(~(I2 == I3));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(I == I)

% compare results
res = true;

% ------------------------------ END OF CODE ------------------------------
