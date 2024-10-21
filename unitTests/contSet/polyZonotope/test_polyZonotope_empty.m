function res = test_polyZonotope_empty
% test_polyZonotope_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_polyZonotope_empty
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

% Authors:       Tobias Ladner
% Written:       16-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D
n = 1;
pZ = polyZonotope.empty(n);
assert(representsa(pZ,'emptySet') && dim(pZ) == 1);

% 5D
n = 5;
pZ = polyZonotope.empty(n);
assert(representsa(pZ,'emptySet') && dim(pZ) == 5);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
