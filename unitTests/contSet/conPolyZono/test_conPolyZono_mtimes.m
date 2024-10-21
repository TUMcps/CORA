function res = test_conPolyZono_mtimes
% test_conPolyZono_mtimes - unit test function of mtimes
%
% Syntax:
%    res = test_conPolyZono_mtimes
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
c = [0;0];
G = [2 0 2; 0 2 2];
GI = [1;6];
E = [1 0 1; 0 1 1; 0 0 0];
A = [2 2 4 -4];
b = 0;
EC = [1 0 1 0; 0 1 1 0; 0 0 0 1];

% multiply with matrix
cPZ = conPolyZono(c,G,E,A,b,EC,GI);
M = [1 2; 3 4];
McPZ = M * cPZ;
assert(compareMatrices(McPZ.c, M*c));
assert(compareMatrices(McPZ.G, M*G));
assert(compareMatrices(McPZ.GI, M*GI));
assert(compareMatrices(McPZ.E, E));
assert(compareMatrices(McPZ.A, A));
assert(compareMatrices(McPZ.b, b));
assert(compareMatrices(McPZ.EC, EC));

% test without GI
cPZ = conPolyZono(c,G,E,A,b,EC);
M = [1 2; 3 4];
McPZ = M * cPZ;

assert(compareMatrices(McPZ.c, M*c));
assert(compareMatrices(McPZ.G, M*G));
assert(isempty(McPZ.GI));

% test with only c
cPZ = conPolyZono(c);
M = [1 2; 3 4];
McPZ = M * cPZ;

assert(compareMatrices(McPZ.c, M*c));
assert(isempty(McPZ.G));
assert(isempty(McPZ.GI));


% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
