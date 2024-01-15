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

res = true(0);

% 2D, empty
C1 = capsule.empty(2);
res(end+1,1) = isemptyobject(C1);

% 2D, bounded
C2 = capsule([1;1],[0;1],0.5);
res(end+1,1) = ~isemptyobject(C2);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
