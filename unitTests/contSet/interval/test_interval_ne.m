function res = test_interval_ne
% test_interval_ne - unit test function of ne, overloaded '~=' operator
%
% Syntax:
%    res = test_interval_ne
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% TEST 1: Analytical ------------------------------------------------------
% create interval
lower1 = [-2; -4; -3];
upper1 = [ 2;  3;  1];
Int1 = interval(lower1, upper1);
lower2 = [-3;  0; -4];
upper2 = [ 2;  3;  1];
Int2 = interval(lower2, upper2);

% compute non-equality
assert(Int1 ~= Int2);    % different intervals
assert(~(Int1 ~= Int1)); % same interval
% -------------------------------------------------------------------------

% TEST 2: Random ----------------------------------------------------------
% create interval
dim = floor(1 + 9*rand(1));
lower1 = -10*rand(dim,1);
upper1 = 10*rand(dim,1);
Int1 = interval(lower1, upper1);
lower2 = -10*rand(dim,1);
upper2 = 10*rand(dim,1);
Int2 = interval(lower2, upper2);

% compute non-equality
res_rand = Int1 ~= Int2;
res_rand_true = any(lower1 ~= lower2) || any(upper1 ~= upper2);
assert(res_rand == res_rand_true);
% -------------------------------------------------------------------------

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(~(I ~= I));

% compare results
res = true;

% ------------------------------ END OF CODE ------------------------------
