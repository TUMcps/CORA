function res = test_conPolyZono_uminus
% test_conPolyZono_uminus - unit test function of uminus
%
% Syntax:
%    res = test_conPolyZono_uminus
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
cPZ = conPolyZono(c,G,E,A,b,EC,GI);

% negate
ncPZ = -cPZ;
assert(all(ncPZ.c == -c, 'all'));
assert(all(ncPZ.G == -G, 'all'));
assert(all(ncPZ.GI == -GI, 'all'));
assert(all(ncPZ.E == E, 'all'));
assert(all(ncPZ.A == A, 'all'));
assert(all(ncPZ.b == b, 'all'));
assert(all(ncPZ.EC == EC, 'all'));

% compare with -1 * cPZ
assert(isequal(ncPZ, -1*cPZ));

% test empty case
assert(isemptyobject(-conPolyZono.empty(2)));

% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
