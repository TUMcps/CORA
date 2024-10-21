function res = test_interval_log
% test_interval_log - unit test function of natural logarithm for intervals
%    overloaded 'log()' function for intervals
%
% Syntax:
%    res = test_interval_log
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
% Written:       07-February-2016
% Last update:   08-June-2020 (MW, rewrite based on new NaN/Inf handling)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-9;

% bounded input, unbounded output
I = interval([0, 1],[2, exp(1)]);
I_log = log(I);
I_true = interval([-Inf,0],[0.6931471805599,1]);
assert(isequal(I_log,I_true,tol));

% out of bounds: throws an error since we cannot instantiate NaN intervals
I = interval([-2; -2],[-1; 0]);
assertThrowsAs(@log,'CORA:outOfDomain',I);

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = abs(interval(lb,ub));
I_log = log(I);
assert(isequal(I_log,interval(log(I.inf),log(I.sup))))


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
