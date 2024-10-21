function res = test_polyZonotope_uplus
% test_polyZonotope_uplus - unit test function of uminus
%
% Syntax:
%    res = test_polyZonotope_uplus
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

% Authors:       Tobias Ladner
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
c = [1;2];
G = [1 0 1; 0 1 1];
GI = [1;2];
E = [1 0 4; 2 1 3];
pZ = polyZonotope(c,G,GI,E);

% plus
ppZ = +pZ;
assert(all(ppZ.c == c, 'all'));
assert(all(ppZ.G == G, 'all'));
assert(all(ppZ.GI == GI, 'all'));

% compare with pZ
assert(isequal(ppZ, pZ));

% test empty case
assert(isemptyobject(+polyZonotope.empty(2)));

% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
