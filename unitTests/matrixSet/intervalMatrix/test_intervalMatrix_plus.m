function res = test_intervalMatrix_plus
% test_intervalMatrix_plus - unit test function for adding an
% interval matrix according to [1]
% 
% 
% Syntax:
%    res = test_intervalMatrix_plus
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
% See also: -

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

%% addition of interval matrices
C_true = intervalMatrix([-1, -1.5; 1.75, 0.5],[3, 4.5; 1.25, 4.5]);
C = A + B;

% initialize partial results
resPartial = [];

% same center?
resPartial(end+1) = (max(max(abs(center(C) - center(C_true)))) == 0);
% same radius?
resPartial(end+1) = (max(max(abs(rad(C) - rad(C_true)))) == 0);
% same infumum?
resPartial(end+1) = (max(max(abs(infimum(C.int) - infimum(C_true.int)))) == 0);
% same supremum?
resPartial(end+1) = (max(max(abs(supremum(C.int) - supremum(C_true.int)))) == 0);

%% addition with a normal matrix
M = [1, 2; -3 1];
C_true = intervalMatrix([1, 2.5; -2.25, -0.5],[1, 1.5; 0.25, 0.5]);
C = M + A;

% same center?
resPartial(end+1) = (max(max(abs(center(C) - center(C_true)))) == 0);
% same radius?
resPartial(end+1) = (max(max(abs(rad(C) - rad(C_true)))) == 0);
% same infumum?
resPartial(end+1) = (max(max(abs(infimum(C.int) - infimum(C_true.int)))) == 0);
% same supremum?
resPartial(end+1) = (max(max(abs(supremum(C.int) - supremum(C_true.int)))) == 0);

%% multiplication with a scalar
m = 3;
C_true = intervalMatrix([3, 3.5; 3.75, 1.5],[1, 1.5; 0.25, 0.5]);
C = m + A;

% same center?
resPartial(end+1) = (max(max(abs(center(C) - center(C_true)))) == 0);
% same radius?
resPartial(end+1) = (max(max(abs(rad(C) - rad(C_true)))) == 0);
% same infumum?
resPartial(end+1) = (max(max(abs(infimum(C.int) - infimum(C_true.int)))) == 0);
% same supremum?
resPartial(end+1) = (max(max(abs(supremum(C.int) - supremum(C_true.int)))) == 0);

%result of all tests
res = all(resPartial);

% ------------------------------ END OF CODE ------------------------------
