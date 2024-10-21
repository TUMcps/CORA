function res = test_interval_reshape
% test_interval_reshape - unit test function of reshape
%
% Syntax:
%    res = test_interval_reshape
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
% Written:       29-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create interval
initSize = [6, 4];
lower = -rand(initSize);
upper = rand(initSize);
I = interval(lower, upper);

% compute reshape
reshapeSize = [8, 3];
Int_reshape = reshape(I, reshapeSize);

% true solution
lower_true = reshape(lower, reshapeSize);
upper_true = reshape(upper, reshapeSize);
Int_true = interval(lower_true, upper_true);

% compare results
assert(Int_reshape == Int_true);

% wrong reshape dimensions should throw error
reshapeSize = [9,2];
assertThrowsAs(@reshape,'MATLAB:getReshapeDims:notSameNumel',I,reshapeSize);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
