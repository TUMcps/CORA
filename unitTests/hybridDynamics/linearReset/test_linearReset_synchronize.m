function res = test_linearReset_synchronize
% test_linearReset_synchronize - test function for the synchronization of
%    a list of linear reset functions
%
% Syntax:
%    res = test_linearReset_synchronize
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
% See also: test_nonlinearReset_synchronize

% Authors:       Mark Wetzlinger
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = eps;

% init linear reset functions
A1 = blkdiag([1 2; 0 -1],zeros(3));
B1 = blkdiag([2 0 1; -1 0 0],zeros(3,1));
c1 = [1; -5; zeros(3,1)];
linReset1 = linearReset(A1,B1,c1);
A2 = blkdiag(zeros(2),[-1 0; 5 3],0);
B2 = blkdiag(zeros(2,3),[-1;0;0]);
c2 = [zeros(2,1); 1; -2; 0];
linReset2 = linearReset(A2,B2,c2);

% synchronize
linResets = [linReset1;linReset2];
idStates = [false;false;false;false;true];
linReset_sync = synchronize(linResets,idStates);

A_true = A1 + A2 + diag([0;0;0;0;1]);
B_true = B1 + B2;
c_true = c1 + c2;
linReset_sync_true = linearReset(A_true,B_true,c_true);

assert(isequal(linReset_sync,linReset_sync_true,tol));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
