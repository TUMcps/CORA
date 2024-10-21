function res = test_conZonotope_uplus
% test_conZonotope_uplus - unit test function of uplus
%
% Syntax:
%    res = test_conZonotope_uplus
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
A = [2 -4 2];
b = 0;
cZ = conZonotope(c,G,A,b);

% plus
pcZ = +cZ;
assert(all(pcZ.c == c, 'all'));
assert(all(pcZ.G == G, 'all'));
assert(all(pcZ.A == A, 'all'));
assert(all(pcZ.b == b, 'all'));

% compare with cZ
assert(isequal(pcZ, cZ));

% test empty case
assert(isemptyobject(+conZonotope.empty(2)));

% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
