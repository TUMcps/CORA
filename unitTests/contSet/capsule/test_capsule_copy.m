function res = test_capsule_copy
% test_capsule_copy - unit test function of copy
%
% Syntax:
%    res = test_capsule_copy
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
% Written:       02-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D capsule
C = capsule([1;-1],[3;1],0.5);
C_copy = copy(C);
assert(isequal(C,C_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
