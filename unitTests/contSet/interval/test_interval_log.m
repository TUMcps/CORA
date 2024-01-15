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

tol = 1e-9;
res = true(0);

% bounded input, unbounded output
I = interval([0, 1],[2, exp(1)]);
I_log = log(I);
I_true = interval([-Inf,0],[0.6931471805599,1]);
res(end+1,1) = isequal(I_log,I_true,tol);

% out of bounds: throws an error since we cannot instantiate NaN intervals
I = interval([-2; -2],[-1; 0]);
try
    log(I);
    res = false;
end


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
