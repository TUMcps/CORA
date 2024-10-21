function res = test_nonlinearReset_nonlinearReset
% test_nonlinearReset_nonlinearReset - test function for nonlinearReset
%    constructor
%
% Syntax:
%    res = test_nonlinearReset_nonlinearReset
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
nonlinReset = nonlinearReset();
assert(nonlinReset.preStateDim == 0);
assert(nonlinReset.inputDim == 1);
assert(nonlinReset.postStateDim == 0);

% only states
f = @(x,u) [x(1)*x(2); x(2)];
nonlinReset = nonlinearReset(f);
assert(nonlinReset.preStateDim == 2);
assert(nonlinReset.inputDim == 1);
assert(nonlinReset.postStateDim == 2);
% ...include dimensions
nonlinReset = nonlinearReset(f,3,2,2);
assert(nonlinReset.preStateDim == 3);
assert(nonlinReset.inputDim == 2);
assert(nonlinReset.postStateDim == 2);
% wrong pre-state dimension, post-state dimension
assertThrowsAs(@nonlinearReset,'CORA:wrongInputInConstructor',f,1,1,2);
assertThrowsAs(@nonlinearReset,'CORA:wrongInputInConstructor',f,2,1,1);

% states and inputs
f = @(x,u) [x(1) - u(1); x(1)*x(2)];
nonlinReset = nonlinearReset(f);
assert(nonlinReset.preStateDim == 2);
assert(nonlinReset.inputDim == 1);
assert(nonlinReset.postStateDim == 2);
% ...include dimensions
nonlinReset = nonlinearReset(f,3,2,2);
assert(nonlinReset.preStateDim == 3);
assert(nonlinReset.inputDim == 2);
assert(nonlinReset.postStateDim == 2);
% wrong pre-state dimension, input dimension, post-state dimension
assertThrowsAs(@nonlinearReset,'CORA:wrongInputInConstructor',f,1,1,2);
assertThrowsAs(@nonlinearReset,'CORA:wrongInputInConstructor',f,2,0,2);
assertThrowsAs(@nonlinearReset,'CORA:wrongInputInConstructor',f,2,1,1);

% states and inputs, different output dimension
f = @(x,u) x(1)*x(2) - u(1);
nonlinReset = nonlinearReset(f);
assert(nonlinReset.preStateDim == 2);
assert(nonlinReset.inputDim == 1);
assert(nonlinReset.postStateDim == 1);
% ...include dimensions
nonlinReset = nonlinearReset(f,3,2,1);
assert(nonlinReset.preStateDim == 3);
assert(nonlinReset.inputDim == 2);
assert(nonlinReset.postStateDim == 1);
% wrong pre-state dimension, input dimension, post-state dimension
assertThrowsAs(@nonlinearReset,'CORA:wrongInputInConstructor',f,1,1,1);
assertThrowsAs(@nonlinearReset,'CORA:wrongInputInConstructor',f,2,0,1);
assertThrowsAs(@nonlinearReset,'CORA:wrongInputInConstructor',f,2,1,0);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
