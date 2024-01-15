function res = test_interval_horzcat
% test_interval_horzcat - unit test function of the operator for
%    horizontal concatenation, e.g. a = [b,c,d];
%
% Syntax:
%    res = test_interval_horzcat
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
% Written:       16-January-2016
% Last update:   04-December-2023 (MW, add unbounded case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true(0);

% bounded
I1 = interval([-5; -4; -3; 0; 0; 5], [-2; 0; 2; 0; 5; 8]);
I2 = I1 + 1;
I3 = I2 + 2;
I = [I1, I2, I3];
I_true = interval([-5 -4 -2; -4 -3 -1; -3 -2 0; 0 1 3; 0 1 3; 5 6 8],...
    [-2 -1 1; 0 1 3; 2 3 5; 0 1 3; 5 6 8; 8 9 11]);
res(end+1,1) = isequal(I,I_true,tol);

% unbounded
I1 = interval([-Inf;-2],[2;3]);
I2 = interval([-2 3; 2 5],[1 Inf; 4 8]);
I = [I1, I2];
I_true = interval([-Inf -2 3; -2 2 5],[2 1 Inf; 3 4 8]);
res(end+1,1) = isequal(I,I_true,tol);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
