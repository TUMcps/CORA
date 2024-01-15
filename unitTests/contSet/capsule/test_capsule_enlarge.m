function res = test_capsule_enlarge
% test_capsule_enlarge - unit test function of enlarge
%
% Syntax:
%    res = test_capsule_enlarge
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
% Written:       15-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate capsule
C = capsule([1;1],[1;1],0.5);
factor = 2;

% compute enlarged capsule
C_enlarged = enlarge(C,factor);

% true solution
C_true = capsule([1;1],[2;2],1);

% compare results
res = isequal(C_true,C_enlarged);

% ------------------------------ END OF CODE ------------------------------
