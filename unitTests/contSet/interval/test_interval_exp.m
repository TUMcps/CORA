function res = test_interval_exp
% test_interval_exp - unit test function of exponential function
%
% Syntax:
%    res = test_interval_exp
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
% See also: mtimes

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       14-January-2016
% Last update:   04-December-2023 (MW, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define problem
tol = 1e-9;
res = true(0);

% bounded
I = interval([-5, -4, -3, 0, 0, 5], [-2, 0, 2, 0, 5, 8]);
I_exp = exp(I);
I_true = interval([0.006737947, 0.0183156389, 0.0497870684, 1, 1, 148.4131591026],...
    [0.135335283, 1, 7.3890560989, 1, 148.4131591026, 2980.9579870418]);
res(end+1,1) = isequal(I_exp,I_true,tol);

% unbounded
I = interval(-Inf,0);
I_exp = exp(I);
I_true = interval(0,1);
res(end+1,1) = isequal(I_exp,I_true,tol);

% unbounded
I = interval(1,Inf);
I_exp = exp(I);
I_true = interval(exp(1),Inf);
res(end+1,1) = isequal(I_exp,I_true,tol);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
