function res = test_interval_transpose
% test_interval_transpose - unit_test_function of transpose,
%    overloaded '.'' operator for intervals
%
% Syntax:
%    res = test_interval_transpose
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
% Last update:   05-December-2023 (MW, add empty and unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true(0);

% empty
I = interval.empty(2);
I_transpose = I.';
res(end+1,1) = representsa(I_transpose,'emptySet') && dim(I) == 2;

% vector
I = interval([-5, -4, -3, 0, 0, 5], [-2, 0, 2, 0, 5, 8]);
I_transpose = I.';
I_true = interval([-5; -4; -3; 0; 0; 5], [-2; 0; 2; 0; 5; 8]);
res(end+1,1) = isequal(I_transpose,I_true,tol);

I = interval([-5; -4; -3; 0; 0; 5], [-2; 0; 2; 0; 5; 8]);
I_transpose = I.';
I_true = interval([-5, -4, -3, 0, 0, 5], [-2, 0, 2, 0, 5, 8]);
res(end+1,1) = isequal(I_transpose,I_true,tol);

% unbounded
I = interval([-Inf; 2; 1],[Inf; Inf; 2]);
I_transpose = I.';
I_true = interval([-Inf, 2, 1],[Inf, Inf, 2]);
res(end+1,1) = isequal(I_transpose,I_true,tol);

% matrix
I = interval([-2 -1; 4 2; -1 0],[1 2; 8 4; 1 0]);
I_transpose = I.';
I_true = interval([-2 4 -1; -1 2 0],[1 8 1; 2 4 0]);
res(end+1,1) = isequal(I_transpose,I_true,tol);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
