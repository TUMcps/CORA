function res = test_capsule_isemptyobject
% test_capsule_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_capsule_isemptyobject
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
% Written:       03-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate capsules
C1 = capsule();
C2 = capsule([1;1],[],0.5);

% check results
res = isemptyobject(C1) && ~isemptyobject(C2);

% ------------------------------ END OF CODE ------------------------------
