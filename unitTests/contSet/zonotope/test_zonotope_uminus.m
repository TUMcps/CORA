function res = test_zonotope_uminus
% test_zonotope_uminus - unit test function of uminus
%
% Syntax:
%    res = test_zonotope_uminus
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
Z = zonotope(c, G);

% negate
nZ = -Z;
assert(all([nZ.c,nZ.G] == -[c,G], 'all'));

% compare with -1 * Z
assert(isequal(nZ, -1*Z));

% test empty case
assert(isemptyobject(-zonotope.empty(2)));

% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
