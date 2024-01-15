function res = test_intervalMatrix_supremum
% test_intervalMatrix_supremum - unit test function for supremum
% 
% Syntax:
%    res = test_intervalMatrix_supremum
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

% Authors:       Mark Wetzlinger
% Written:       16-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = eps;

% 2x2 interval matrix
C = [2 3; 1 2]; D = [1 0; 1 1];
intMat = intervalMatrix(C,D);
lb = supremum(intMat);
lb_true = C+D;
res(end+1,1) = all(all(withinTol(lb,lb_true,tol)));

% 1x3 interval matrix
C = [2; -1; 1]; D = [1; 2; 1];
intMat = intervalMatrix(C,D);
lb = supremum(intMat);
lb_true = C+D;
res(end+1,1) = all(all(withinTol(lb,lb_true,tol)));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
