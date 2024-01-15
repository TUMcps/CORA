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

% empty polynomial zonotope
pZ = polyZonotope.empty(2);
res(end+1,1) = isemptyobject(pZ);

% instantiate polynomial zonotopes
pZ = polyZonotope([0;0],[0 4 1 -1 2; 1 2 -1 -1 1],...
    [-7 1 1;15 1 -1],[1 0 0 0 1;0 1 0 3 2; 0 0 1 1 0]);
res(end+1,1) = ~isemptyobject(pZ);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
