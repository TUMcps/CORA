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
res_e = or(Int1,I_empty) == Int1;
% compare results
res_analytical = all(infimum(IntUnion) == infimum(IntUnion_true)) && ...
        all(supremum(IntUnion) == supremum(IntUnion_true)) && res_e;
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
res_rand = all(infimum(IntUnion) == infimum(IntUnion_true)) && ...
        all(supremum(IntUnion) == supremum(IntUnion_true));
% -------------------------------------------------------------------------

% add results
res = res_analytical && res_rand;

% ------------------------------ END OF CODE ------------------------------
