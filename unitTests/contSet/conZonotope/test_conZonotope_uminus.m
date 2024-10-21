function res = test_conZonotope_uminus
% test_conZonotope_uminus - unit test function of uminus
%
% Syntax:
%    res = test_conZonotope_uminus
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

% negate
ncZ = -cZ;
assert(all(ncZ.c == -c, 'all'));
assert(all(ncZ.G == -G, 'all'));
assert(all(ncZ.A == A, 'all'));
assert(all(ncZ.b == b, 'all'));

% compare with -1 * cZ
assert(isequal(ncZ, -1*cZ));

% test empty case
assert(isemptyobject(-conZonotope.empty(2)));

% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
