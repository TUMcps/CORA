function res = test_polyZonotope_cartProd
% test_polyZonotope_cartProd - unit test function for the Cartesian product
%    of a polynomial zonotope and another set or point
%
% Syntax:
%    res = test_polyZonotope_cartProd
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
% Written:       04-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% different polynomial zonotopes
c = [1; 2];
G = [1 2 1 -3; 1 -1 2 -1];
G_def = zeros(2,4);
E = [1 0 0 2; 0 1 2 1];
E_def = eye(4);
GI = [1; 0];
GI_def = zeros(2,1);
id = [5;6];
id_def2 = [1;2];
id_def4 = [1;2;3;4];

pZ_c = polyZonotope(c);
pZ_cG = polyZonotope(c,G);
pZ_cGE = polyZonotope(c,G,[],E);
pZ_cGEid = polyZonotope(c,G,[],E,id);
pZ_all = polyZonotope(c,G,GI,E,id);


%%% polyZonotope x polyZonotope

% center x center
pZ = cartProd(pZ_c,pZ_c);
assert(all(pZ.c == [c;c]));
assert(isempty(pZ.G) && isempty(pZ.GI) && isempty(pZ.E) && isempty(pZ.id));

% center x center, dependent generators
pZ = cartProd(pZ_c,pZ_cG);
assert(all(pZ.c == [c;c]));
assert(all(pZ.G == [G_def;G],'all'));
assert(isempty(pZ.GI));
assert(all(pZ.E == E_def,'all'));
assert(all(pZ.id == id_def4));
% vice versa
pZ = cartProd(pZ_cG,pZ_c);
assert(all(pZ.c == [c;c]));
assert(all(pZ.G == [G;G_def],'all'));
assert(isempty(pZ.GI));
assert(all(pZ.E == E_def,'all'));
assert(all(pZ.id == id_def4));

% center x center, dependent generators, exponent matrix
pZ = cartProd(pZ_c,pZ_cGE);
assert(all(pZ.c == [c;c]));
assert(all(pZ.G == [G_def;G],'all'));
assert(isempty(pZ.GI));
assert(all(pZ.E == E,'all'));
assert(all(pZ.id == id_def2));
% vice versa
pZ = cartProd(pZ_cGE,pZ_c);
assert(all(pZ.c == [c;c]));
assert(all(pZ.G == [G;G_def],'all'));
assert(isempty(pZ.GI));
assert(all(pZ.E == E,'all'));
assert(all(pZ.id == id_def2));

% center x center, dependent generators, exponent matrix, identifier vector
pZ = cartProd(pZ_c,pZ_cGEid);
assert(all(pZ.c == [c;c]));
assert(all(pZ.G == [G_def;G],'all'));
assert(isempty(pZ.GI));
assert(all(pZ.E == E,'all'));
assert(all(pZ.id == id));
% vice versa
pZ = cartProd(pZ_cGEid,pZ_c);
assert(all(pZ.c == [c;c]));
assert(all(pZ.G == [G;G_def],'all'));
assert(isempty(pZ.GI));
assert(all(pZ.E == E,'all'));
assert(all(pZ.id == id));

% center x all
pZ = cartProd(pZ_c,pZ_all);
assert(all(pZ.c == [c;c]));
assert(all(pZ.G == [G_def;G],'all'));
assert(all(pZ.GI == [GI_def;GI]));
assert(all(pZ.E == E,'all'));
assert(all(pZ.id == id));
% vice versa
pZ = cartProd(pZ_all,pZ_c);
assert(all(pZ.c == [c;c]));
assert(all(pZ.G == [G;G_def],'all'));
assert(all(pZ.GI == [GI;GI_def]));
assert(all(pZ.E == E,'all'));
assert(all(pZ.id == id));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
