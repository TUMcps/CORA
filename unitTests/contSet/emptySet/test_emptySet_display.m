function res = test_emptySet_display
% test_emptySet_display - unit test function of display
%
% Syntax:
%    res = test_emptySet_display
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

% 0-dimensional empty set
O = emptySet(0)

% n-dimensional empty set
n = 2;
O = emptySet(n)

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
