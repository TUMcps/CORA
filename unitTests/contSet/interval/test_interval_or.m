function res = test_interval_or
% test_interval_or - unit test function of or
%
% Syntax:
%    res = test_interval_or
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
% Last update:   16-September-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% TEST 1: Analytical ------------------------------------------------------
% create intervals
lower1 = [-2; -4; -3];
upper1 = [ 2;  3;  1];
lower2 = [-1; -6; -3];
upper2 = [ 1;  5;  2];
Int1 = interval(lower1, upper1);
Int2 = interval(lower2, upper2);

% compute union
IntUnion = or(Int1, Int2);

% true size
lower_true = [-2; -6; -3];
upper_true = [ 2;  5;  2];
IntUnion_true = interval(lower_true, upper_true);

% empty set
I_empty = interval.empty(3);
assert(or(Int1,I_empty) == Int1);
% compare results
assert(all(infimum(IntUnion) == infimum(IntUnion_true)))
assert(all(supremum(IntUnion) == supremum(IntUnion_true)))

% -------------------------------------------------------------------------

% TEST 2: Random ----------------------------------------------------------
% create random intervals
dim = floor(1 + 9*rand(1));
lower1 = -10*rand(dim,1);
upper1 = 10*rand(dim,1);
lower2 = -10*rand(dim,1);
upper2 = 10*rand(dim,1);
Int1 = interval(lower1, upper1);
Int2 = interval(lower2, upper2);

% compute union
IntUnion = or(Int1, Int2);

% true size
lower_true = min(lower1, lower2);
upper_true = max(upper1, upper2);
IntUnion_true = interval(lower_true, upper_true);

% compare results
assert(all(infimum(IntUnion) == infimum(IntUnion_true)))
assert(all(supremum(IntUnion) == supremum(IntUnion_true)))

% TEST 3: Empty sets ------------------------------------------------------

I = interval([1;2],[3;4]);
Iempty = interval.empty(2);

I_or = I | Iempty;
assert(isequal(I,I_or));

I_or = Iempty | I;
assert(isequal(I,I_or));

% -------------------------------------------------------------------------

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I1 = interval(lb,ub);
I2 = I1+1;
I = I1 | I2;
I_true = interval(lb,ub+1);

assert(isequal(I,I_true))

% -------------------------------------------------------------------------

% check if broadcasting works
w = [1 2 3; 4 5 6];
W = [0;0];

% check all combinations
I = interval(w) | W;
assert(contains(I,w) && contains(I,W))
I = W | interval(w);
assert(contains(I,w) && contains(I,W))
I = interval(W) | w;
assert(contains(I,w) && contains(I,W))
I = w | interval(W);
assert(contains(I,w) && contains(I,W))

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
