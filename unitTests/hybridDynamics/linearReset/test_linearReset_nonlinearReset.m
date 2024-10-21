function res = test_linearReset_nonlinearReset
% test_linearReset_nonlinearReset - test function for conversion of a
%    linearReset object into a nonlinearReset object
%
% Syntax:
%    res = test_linearReset_nonlinearReset
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
% Written:       08-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty
linReset = linearReset();
nonlinReset = nonlinearReset(linReset);
assert(nonlinReset.preStateDim == 0);
assert(nonlinReset.inputDim == 1);
assert(nonlinReset.postStateDim == 0);

% only state matrix
A = [1 0; 0 1];
linReset = linearReset(A);
nonlinReset = nonlinearReset(linReset);
assert(nonlinReset.preStateDim == 2);
assert(nonlinReset.inputDim == 1);
assert(nonlinReset.postStateDim == 2);
f = @(x,u) [x(1); x(2)];
nonlinReset_ = nonlinearReset(f);
assert(isequal(nonlinReset,nonlinReset_));

% state matrix 1D
A = 2;
linReset = linearReset(A);
nonlinReset = nonlinearReset(linReset);
assert(nonlinReset.preStateDim == 1);
assert(nonlinReset.inputDim == 1);
assert(nonlinReset.postStateDim == 1);
f = @(x,u) 2*x(1);
nonlinReset_ = nonlinearReset(f);
assert(isequal(nonlinReset,nonlinReset_));

% state matrix and input matrix
A = [1 0; 0 1]; B = [1; -1];
linReset = linearReset(A,B);
nonlinReset = nonlinearReset(linReset);
assert(nonlinReset.preStateDim == 2);
assert(nonlinReset.inputDim == 1);
assert(nonlinReset.postStateDim == 2);
f = @(x,u) [x(1) + u(1); x(2) - u(1)];
nonlinReset_ = nonlinearReset(f);
assert(isequal(nonlinReset,nonlinReset_));

% state matrix and offset
A = [1 0; 0 1]; c = [1; -1];
linReset = linearReset(A,[],c);
nonlinReset = nonlinearReset(linReset);
assert(nonlinReset.preStateDim == 2);
assert(nonlinReset.inputDim == 1);
assert(nonlinReset.postStateDim == 2);
f = @(x,u) [x(1) + 1; x(2) - 1];
nonlinReset_ = nonlinearReset(f);
assert(isequal(nonlinReset,nonlinReset_));

% state matrix, input matrix, and offset
A = [1 0; 0 1]; B = [0;1]; c = [1; -1];
linReset = linearReset(A,B,c);
nonlinReset = nonlinearReset(linReset);
assert(nonlinReset.preStateDim == 2);
assert(nonlinReset.inputDim == 1);
assert(nonlinReset.postStateDim == 2);
f = @(x,u) [x(1) + 1; x(2) + u(1) - 1];
nonlinReset_ = nonlinearReset(f);
assert(isequal(nonlinReset,nonlinReset_));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
