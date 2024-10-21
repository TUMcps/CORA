function res = test_emptySet_lift
% test_emptySet_lift - unit test function of lift
%
% Syntax:
%    res = test_emptySet_lift
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init fullspace
n = 4;
O = emptySet(n);

% higher-dimensional space
N = 6;
dims = [1,2,5,6];
O_ = lift(O,N,dims);
fs_true = emptySet(N);
assert(isequal(fs_true,O_));

% not enough specified dimensions: fix error message before this...
% projDims = [1,2,5];
% try
%     O_ = lift(O,N,projDims);
%     res(end+1,1) = false;
% catch
%     res(end+1,1) = true;
% end

% dimensions out of range
projDims = [-1,2,3,5];
assertThrowsAs(@lift,'CORA:wrongValue',O,N,projDims);

% higher-dimensional space smaller than original space
N = 3;
projDims = [1,2,5,6];
assertThrowsAs(@lift,'CORA:wrongValue',O,N,projDims);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
