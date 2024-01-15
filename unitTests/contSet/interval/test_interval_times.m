function res = test_interval_times
% test_interval_times - unit test function of times,
%    overloaded '.*' operator for intervals
%
% Syntax:
%    res = test_interval_times
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
%                05-December-2023 (MW, add unbounded and matrix cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true(0);

% numeric, interval
I = interval([2;-2],[4; 5]);
num = -2;
I_times = num .* I;
I_true = interval([-8;-10],[-4;4]);
res(end+1,1) = isequal(I_times,I_true,tol);
I_times = I .* num;
res(end+1,1) = isequal(I_times,I_true,tol);

% unbonded interval, numeric
I = interval([-Inf;-2],[4; Inf]);
num = -2;
I_times = I .* num;
I_true = interval([-8;-Inf],[Inf;4]);
res(end+1,1) = isequal(I_times,I_true,tol);

% interval matrix, numeric
I = interval([-2 -1; 4 3],[5 3; 9 8]);
num = 4;
I_times = I .* num;
I_true = interval([-8 -4; 16 12],[20 12; 36 32]);
res(end+1,1) = isequal(I_times,I_true,tol);
I_times = num .* I;
res(end+1,1) = isequal(I_times,I_true,tol);

% interval, interval
I1 = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
I2 = interval([-6.1, -4.5, -3.3, 0, 0, 5], [-2.2, 0.0, 2.8, 0, 5.7, 8.2]);
I_times = I1 .* I2;
I_true = interval([4.4,0,-8.4,0,0,25],[30.5,18,9.9,0,28.5,65.6]);
res(end+1,1) = isequal(I_times,I_true,tol);

% interval matrix, interval
I1 = interval([-2 1; 1 0; -4 -3],[3 2; 6 4; -2 0]);
I2 = interval([-1; 0; 1], [1; 1; 2]);
I_times = I1 .* I2;
I_true = interval([-3 -2; 0 0; -8 -6],[3 2; 6 4; -2 0]);
res(end+1,1) = isequal(I_times,I_true,tol);
I_times = I2 .* I1;
res(end+1,1) = isequal(I_times,I_true,tol);

% interval matrix, interval matrix
I1 = interval([-2 1; 3 2],[0 2; 4 3]);
I2 = interval([-1 -4; 0 1],[1 -3; 1 5]);
I_times = I1 .* I2;
I_true = interval([-2 -8; 0 2],[2 -3; 4 15]);
res(end+1,1) = isequal(I_times,I_true,tol);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
