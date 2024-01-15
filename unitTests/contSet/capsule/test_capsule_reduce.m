function res = test_capsule_reduce
% test_capsule_reduce - unit test function of reduce
%
% Syntax:
%    res = test_capsule_reduce
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty capsule
C = capsule.empty(2);
C_reduce = reduce(C);
res(end+1,1) = isequal(C,C_reduce);

% 2D capsule
C = capsule([1;4],[1;-2],2);
C_ = reduce(C);
res(end+1,1) = isequal(C,C_);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
