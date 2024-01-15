function res = test_interval_mpower
% test_interval_mpower - unit_test_function of power,
%    overloaded '^' operator for intervals
%
% Syntax:
%    res = test_interval_mpower
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

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       05-January-2016
% Last update:   04-December-2023 (MW, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true(0);

% bounded
I = interval(0, 2);
I_mpower = I ^ 1;
res(end+1,1) = isequal(I_mpower,I,tol);

I = interval(0, 2);
I_mpower = I ^ 2;
I_true = interval(0,4);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(-2, 0);
I_mpower = I ^ 2;
I_true = interval(0,4);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(-2, 0);
I_mpower = I ^ 3;
I_true = interval(-8,0);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(-3, 2);
I_mpower = I ^ 2;
I_true = interval(0,9);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(-3, 2);
I_mpower = I ^ 3;
I_true = interval(-27,8);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(-3, -2);
I_mpower = I ^ 2;
I_true = interval(4,9);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(-3, -2);
I_mpower = I ^ 3;
I_true = interval(-27,-8);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(2, 3);
I_mpower = I ^ 2;
I_true = interval(4,9);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(2, 3);
I_mpower = I ^ 3;
I_true = interval(8,27);
res(end+1,1) = isequal(I_mpower,I_true,tol);

% unbounded
I = interval(-Inf,2);
I_mpower = I ^ 1;
res(end+1,1) = isequal(I_mpower,I,tol);

I = interval(-Inf,2);
I_mpower = I ^ 2;
I_true = interval(0,Inf);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(-Inf,2);
I_mpower = I ^ 3;
I_true = interval(-Inf,8);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(2,Inf);
I_mpower = I ^ 1;
res(end+1,1) = isequal(I_mpower,I,tol);

I = interval(2,Inf);
I_mpower = I ^ 2;
I_true = interval(4,Inf);
res(end+1,1) = isequal(I_mpower,I_true,tol);

I = interval(2,Inf);
I_mpower = I ^ 3;
I_true = interval(8,Inf);
res(end+1,1) = isequal(I_mpower,I_true,tol);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
