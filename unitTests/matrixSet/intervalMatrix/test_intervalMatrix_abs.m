function res = test_intervalMatrix_abs
% test_intervalMatrix_abs - unit test function for absloute values of an
% interval matrix according to [1]
% 
% 
% Syntax:
%    res = test_intervalMatrix_abs
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

%% create interval matrix
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

%% compute absolute value
C_true = [1, 2; 1, 2];
C = abs(A);

% same result?
res = (max(max(abs(C - C_true))) == 0);

% ------------------------------ END OF CODE ------------------------------
