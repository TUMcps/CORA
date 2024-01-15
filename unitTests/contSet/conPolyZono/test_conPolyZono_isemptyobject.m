function res = test_conPolyZono_isemptyobject
% test_conPolyZono_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_conPolyZono_isemptyobject
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

% empty constrained polynomial zonotope
cPZ = conPolyZono.empty(2);
res(end+1,1) = isemptyobject(cPZ);

% 2D constrained polynomial zonotope
c = [0;0];
G = [1 0;0 1];
E = [1 0;0 1];
A = [1 -1];
b = 0;
EC = [2 0;0 1];
cPZ = conPolyZono(c,G,E,A,b,EC);
res(end+1,1) = ~isemptyobject(cPZ);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
