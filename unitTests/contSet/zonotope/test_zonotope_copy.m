function res = test_zonotope_copy
% test_zonotope_copy - unit test function of copy
%
% Syntax:
%    res = test_zonotope_copy
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
% Written:       02-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D zonotope
Z = zonotope(2,4);
Z_copy = copy(Z);
assert(isequal(Z,Z_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
