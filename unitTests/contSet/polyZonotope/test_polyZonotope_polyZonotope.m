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
% Last update:   04-October-2024 (MW, check default properties)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test all different syntaxes from constructor
res = true;

% empty polyZonotope
n = 2;
pZ = polyZonotope.empty(n);
assert(all(size(pZ.c) == [n,0]));
assert(all(size(pZ.G) == [n,0]));
assert(all(size(pZ.GI) == [n,0]));
assert(all(size(pZ.E) == [0,0]));
assert(all(size(pZ.id) == [0,1]));
assert(representsa(pZ,'emptySet'));

% create polynomial zonotope
c = [0;0];
G = [2 0 1;0 2 1];
GI = [0;0.5];
E = [1 0 3;0 1 1];
E_def = eye(3);
id = [5;6];
id_def2 = [1;2];
id_def3 = [1;2;3];

% only center
pZ = polyZonotope(c);
assert(all(size(pZ.G) == [2,0]));

% only center and dependent generator matrix
pZ = polyZonotope(c,G);
assert(all(size(pZ.GI) == [2,0]));
assert(all(pZ.E == E_def,'all'));
assert(all(pZ.id == id_def3));

% center and both generator matrices
pZ = polyZonotope(c,G,GI);
assert(all(pZ.E == E_def,'all'));
assert(all(pZ.id == id_def3));

% only independent generator matrix
pZ = polyZonotope(c,[],GI);
assert(all(size(pZ.G) == [2,0]));
assert(all(size(pZ.E) == [0,0]));
assert(all(size(pZ.id) == [0,1]));

% both generator matrices and exponent matrix
pZ = polyZonotope(c,G,GI,E);
assert(all(pZ.id == id_def2));

% no independent generator matrix
pZ = polyZonotope(c,G,[],E);
assert(all(size(pZ.GI) == [2,0]));
assert(all(pZ.id == id_def2));

% no independent generator matrix, with identifiers
pZ = polyZonotope(c,G,[],E,id);
assert(all(size(pZ.GI) == [2,0]));

% all input arguments
pZ = polyZonotope(c,G,GI,E,id);

% copy constructor
pZ = polyZonotope(pZ);

% ------------------------------ END OF CODE ------------------------------
