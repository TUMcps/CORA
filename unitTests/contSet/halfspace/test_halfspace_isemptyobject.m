function res = test_halfspace_isemptyobject
% test_halfspace_isemptyobject - unit test function of isemptyobject
%
% Syntax:  
%    res = test_halfspace_isemptyobject
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      03-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate halfspaces
hs1 = halfspace();

a = [3; 2; -1];
b = 0.5;
hs2 = halfspace(a,b);

% check results
res = isemptyobject(hs1) && ~isemptyobject(hs2);

%------------- END OF CODE --------------
