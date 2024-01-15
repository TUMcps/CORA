function res = test_intervalMatrix_mtimes
% test_intervalMatrix_mtimes - unit test function for multiplying an
% interval matrix according to [1]
% 
% 
% Syntax:
%    res = test_intervalMatrix_mtimes
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] Althoff, M.; Grebenyuk, D.: Implementation of Interval Arithmetic 
%        in CORA 2016. Proc. of the 3rd International Workshop on Applied 
%        Verification for Continuous and Hybrid Systems, 2016, 91-105
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       06-May-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% create first interval matrix
% center
matrixCenter = [ ...
0, 0.5; ...
0.75, -1.5];

% delta
matrixDelta = [ ...
1, 1.5; ...
0.25, 0.5];

% instantiate interval matrix
A = intervalMatrix(matrixCenter, matrixDelta);

%% create second interval matrix
% center
matrixCenter = [ ...
-1, -2; ...
1, 2];

% delta
matrixDelta = [ ...
2, 3; ...
1, 4];

% instantiate interval matrix
B = intervalMatrix(matrixCenter, matrixDelta);

%% multiplication of interval matrices
C_true = intervalMatrix([1, 3; -3, -6],[6, 14; 4, 11]);
C = A*B;

% initialize partial results
res = true(0);

% same center?
res(end+1,1) = (max(max(abs(center(C) - center(C_true)))) == 0);
% same radius?
res(end+1,1) = (max(max(abs(rad(C) - rad(C_true)))) == 0);
% same infumum?
res(end+1,1) = (max(max(abs(infimum(C.int) - infimum(C_true.int)))) == 0);
% same supremum?
res(end+1,1) = (max(max(abs(supremum(C.int) - supremum(C_true.int)))) == 0);

%% multiplication with a normal matrix
M = [1, 2; -3 1];
C_true = intervalMatrix([1.5, -2.5; 0.75, -3],[1.5, 2.5; 3.25, 5]);
C = M*A;

% same center?
res(end+1,1) = (max(max(abs(center(C) - center(C_true)))) == 0);
% same radius?
res(end+1,1) = (max(max(abs(rad(C) - rad(C_true)))) == 0);
% same infumum?
res(end+1,1) = (max(max(abs(infimum(C.int) - infimum(C_true.int)))) == 0);
% same supremum?
res(end+1,1) = (max(max(abs(supremum(C.int) - supremum(C_true.int)))) == 0);

%% multiplication with a scalar
m = 3;
C_true = intervalMatrix([0, 1.5; 2.25, -4.5],[3, 4.5; 0.75, 1.5]);
C = m*A;

% same center?
res(end+1,1) = (max(max(abs(center(C) - center(C_true)))) == 0);
% same radius?
res(end+1,1) = (max(max(abs(rad(C) - rad(C_true)))) == 0);
% same infumum?
res(end+1,1) = (max(max(abs(infimum(C.int) - infimum(C_true.int)))) == 0);
% same supremum?
res(end+1,1) = (max(max(abs(supremum(C.int) - supremum(C_true.int)))) == 0);

%result of all tests
res = all(res);

% ------------------------------ END OF CODE ------------------------------
