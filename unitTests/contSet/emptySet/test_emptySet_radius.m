function res = test_emptySet_radius
% test_emptySet_radius - unit test function of radius
%
% Syntax:
%    res = test_emptySet_radius
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
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% init empty set
n = 2;
O = emptySet(n);

% compute radius
r = radius(O);
assert(r == 0);

% ------------------------------ END OF CODE ------------------------------
