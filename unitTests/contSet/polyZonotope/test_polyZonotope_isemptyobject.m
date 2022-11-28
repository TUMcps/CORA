function res = test_polyZonotope_isemptyobject
% test_polyZonotope_isemptyobject - unit test function of isemptyobject
%
% Syntax:  
%    res = test_polyZonotope_isemptyobject
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

% instantiate polynomial zonotopes
pZ1 = polyZonotope();
pZ2 = polyZonotope([0;0],[0 4 1 -1 2; 1 2 -1 -1 1],...
    [-7 1 1;15 1 -1],[1 0 0 0 1;0 1 0 3 2; 0 0 1 1 0]);

% check results
res = isemptyobject(pZ1) && ~isemptyobject(pZ2);

%------------- END OF CODE --------------
