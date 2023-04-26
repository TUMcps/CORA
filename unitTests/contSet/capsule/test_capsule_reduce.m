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

% Author:       Mark Wetzlinger
% Written:      23-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate capsule
C = capsule([1;4],[1;-2],2);

% reduce (should have no effect)
C_ = reduce(C);

% compare
res = isequal(C,C_);

%------------- END OF CODE --------------