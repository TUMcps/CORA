function res = test_interval_isBounded
% test_interval_isBounded - unit test function of isBounded
%
% Syntax:
%    res = test_interval_isBounded
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

% Authors:       Tobias Ladner
% Written:       26-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate bounded intervals
I = interval.empty(2);
assert(isBounded(I));

I = interval(1);
assert(isBounded(I));

I = interval([-2;-1;-3],[1;1;2]);
assert(isBounded(I));

% instantiate unbounded intervals
I = interval([-2;-Inf;-3],[1;1;2]);
assert(~isBounded(I));

I = interval([-2;-Inf;-3],[1;Inf;2]);
assert(~isBounded(I));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(isBounded(I));

% check results
res = true;

% ------------------------------ END OF CODE ------------------------------
