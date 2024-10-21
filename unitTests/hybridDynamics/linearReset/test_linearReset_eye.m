function res = test_linearReset_eye
% test_linearReset_eye - test function for instantiation of identity reset
%
% Syntax:
%    res = test_linearReset_eye
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
% Written:       15-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty states
n = 0;
linReset = linearReset.eye(n);
assert(linReset.preStateDim == n);
assert(linReset.inputDim == 1);
assert(linReset.postStateDim == n);

% non-empty states
n = 2;
linReset = linearReset.eye(n);
assert(linReset.preStateDim == n);
assert(linReset.inputDim == 1);
assert(linReset.postStateDim == n);

% states and inputs
n = 3; m = 2;
linReset = linearReset.eye(n,m);
assert(linReset.preStateDim == n);
assert(linReset.inputDim == m);
assert(linReset.postStateDim == n);

% wrong calls
assertThrowsAs(@linearReset.eye,'MATLAB:narginchk:notEnoughInputs');
assertThrowsAs(@linearReset.eye,'MATLAB:narginchk:tooManyInputs',2,2,2);
assertThrowsAs(@linearReset.eye,'CORA:wrongValue',-1);
assertThrowsAs(@linearReset.eye,'CORA:wrongValue',1.5);
assertThrowsAs(@linearReset.eye,'CORA:wrongValue','two');


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
