function res = test_abstractReset_reset
% test_abstractReset_reset - test function for abstractReset constructor
%
% Syntax:
%    res = test_abstractReset_reset
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
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% standard case
n = 1; m = 2; n_ = 1;
reset = abstractReset(n,m,n_);
assert(reset.preStateDim == n);
assert(reset.inputDim == m);
assert(reset.postStateDim == n_);

% copy constructor
n = 1; m = 2; n_ = 1;
reset = abstractReset(n,m,n_);
reset_ = abstractReset(reset);
assert(reset.preStateDim == reset_.preStateDim);
assert(reset.inputDim == reset_.inputDim);
assert(reset.postStateDim == reset_.postStateDim);

% not enough or too many input arguments
assertThrowsAs(@abstractReset,'CORA:numInputArgsConstructor');
assertThrowsAs(@abstractReset,'CORA:numInputArgsConstructor',1,2);
assertThrowsAs(@abstractReset,'CORA:numInputArgsConstructor',1,2,3,4);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
