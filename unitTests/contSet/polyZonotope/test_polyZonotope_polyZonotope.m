function res = test_polyZonotope_polyZonotope
% test_polyZonotope_polyZonotope - unit test function for constructor
%
% Syntax:
%    res = test_polyZonotope_polyZonotope
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
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test all different syntaxes from constructor
res = true;

% empty polyZonotope
pZ = polyZonotope.empty(2);

% create polynomial zonotope
c = [0;0];
G = [2 0 1;0 2 1];
GI = [0;0.5];
E = [1 0 3;0 1 1];
id = [1;2];

% only center and dependent generator matrix
pZ = polyZonotope(c,G);

% copy constructor
pZ = polyZonotope(pZ);

% center and both generator matrices
pZ = polyZonotope(c,G,GI);

% only independent generator matrix
pZ = polyZonotope(c,[],GI);

% both generator matrices and exponent matrix
pZ = polyZonotope(c,G,GI,E);

% no independent generator matrix
pZ = polyZonotope(c,G,[],E);

% all input arguments
pZ = polyZonotope(c,G,GI,E,id);

% no independent generator matrix, with identifiers
pZ = polyZonotope(c,G,[],E,id);

% ------------------------------ END OF CODE ------------------------------
