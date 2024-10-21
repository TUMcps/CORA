function res = test_conPolyZono_empty
% test_conPolyZono_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_conPolyZono_empty
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
% Written:       09-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate constrained polynomial zonotope
n = 2;
cPZ = conPolyZono.empty(n);
assert(isempty(cPZ.c) && isnumeric(cPZ.c) && size(cPZ.c,1) == n);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
