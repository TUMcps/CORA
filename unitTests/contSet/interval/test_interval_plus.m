function res = test_interval_plus
% test_interval_plus - unit_test_function of plus,
%    overloaded '+' operator for intervals
%
% Syntax:
%    res = test_interval_plus
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
% Written:       04-January-2016
% Last update:   13-January-2016 (DG)
%                05-December-2023 (MW, add unbounded, matrix cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true(0);

% empty
I = interval.empty(1);
v = 1;
I_plus = I + v;
res(end+1,1) = representsa(I_plus,'emptySet');

% bounded interval, numeric
v = [2;1];
I = interval([-2;-1],[3; 4]);
I_plus = v + I;
I_true = interval([0;0],[5;5]);
res(end+1,1) = isequal(I_plus,I_true,tol);
I_plus = I + v;
res(end+1,1) = isequal(I_plus,I_true,tol);

% unbounded interval, numeric
I = interval(-Inf,2);
v = 1;
I_plus = I + v;
I_true = interval(-Inf,3);
res(end+1,1) = isequal(I_plus,I_true,tol);
I_plus = v + I;
res(end+1,1) = isequal(I_plus,I_true,tol);

% bounded interval, bounded interval
I1 = interval([-2;-1],[3;4]);
I2 = interval([-1;-3],[1;-1]);
I_plus = I1 + I2;
I_true = interval([-3;-4],[4;3]);
res(end+1,1) = isequal(I_plus,I_true,tol);
I_plus = I2 + I1;
res(end+1,1) = isequal(I_plus,I_true,tol);

% unbounded interval, unbounded interval
I1 = interval([-Inf;-2],[2;4]);
I2 = interval([-1;0],[1;Inf]);
I_plus = I1 + I2;
I_true = interval([-Inf;-2],[3;Inf]);
res(end+1,1) = isequal(I_plus,I_true,tol);
I_plus = I2 + I1;
res(end+1,1) = isequal(I_plus,I_true,tol);

% interval matrix, numeric
I = interval([-2 -1; 0 2],[3 5; 2 3]);
v = 2;
I_plus = I + v;
I_true = interval([0 1; 2 4],[5 7; 4 5]);
res(end+1,1) = isequal(I_plus,I_true,tol);
I_plus = v + I;
res(end+1,1) = isequal(I_plus,I_true,tol);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
