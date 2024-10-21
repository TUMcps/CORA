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

% bounded
I1 = interval([-5; -4; -3; 0; 0; 5], [-2; 0; 2; 0; 5; 8]);
I2 = I1 + 1;
I3 = I2 + 2;
I = [I1, I2, I3];
I_true = interval([-5 -4 -2; -4 -3 -1; -3 -2 0; 0 1 3; 0 1 3; 5 6 8],...
    [-2 -1 1; 0 1 3; 2 3 5; 0 1 3; 5 6 8; 8 9 11]);
assert(isequal(I,I_true,tol));

% unbounded
I1 = interval([-Inf;-2],[2;3]);
I2 = interval([-2 3; 2 5],[1 Inf; 4 8]);
I = [I1, I2];
I_true = interval([-Inf -2 3; -2 2 5],[2 1 Inf; 3 4 8]);
assert(isequal(I,I_true,tol));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
I = [I I];
assert(isequal(size(I),[2,4,2,3]))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
