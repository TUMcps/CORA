function res = test_zonotope_empty
% test_zonotope_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_zonotope_empty
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
Z = zonotope.empty(n);
assert(representsa(Z,'emptySet') && dim(Z) == 1);

% 5D
n = 5;
Z = zonotope.empty(n);
assert(representsa(Z,'emptySet') && dim(Z) == 5);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
