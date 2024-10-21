function res = test_linearReset_linearReset
% test_linearReset_linearReset - test function for linearReset constructor
%
% Syntax:
%    res = test_linearReset_linearReset
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
assert(linReset.preStateDim == 0);
assert(linReset.inputDim == 1);
assert(linReset.postStateDim == 0);

% only state matrix
A = [1 0; 0 1];
linReset = linearReset(A);
assert(all(linReset.A == A,'all'));
assert(all(linReset.B == [0;0]));
assert(all(linReset.c == [0;0]));
assert(linReset.preStateDim == 2);
assert(linReset.inputDim == 1);
assert(linReset.postStateDim == 2);

% state matrix and input matrix
A = [1 0; 0 1]; B = [1; -1];
linReset = linearReset(A,B);
assert(all(linReset.A == A,'all'));
assert(all(linReset.B == B));
assert(all(linReset.c == [0;0]));
assert(linReset.preStateDim == 2);
assert(linReset.inputDim == 1);
assert(linReset.postStateDim == 2);

% state matrix and offset
A = [1 0; 0 1]; c = [1; -1];
linReset = linearReset(A,[],c);
assert(all(linReset.A == A,'all'));
assert(all(linReset.B == [0;0]));
assert(all(linReset.c == c));
assert(linReset.preStateDim == 2);
assert(linReset.inputDim == 1);
assert(linReset.postStateDim == 2);

% state matrix, input matrix, and offset
A = [1 0; 0 1]; B = [0;1]; c = [1; -1];
linReset = linearReset(A,B,c);
assert(all(linReset.A == A,'all'));
assert(all(linReset.B == B));
assert(all(linReset.c == c));
assert(linReset.preStateDim == 2);
assert(linReset.inputDim == 1);
assert(linReset.postStateDim == 2);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
