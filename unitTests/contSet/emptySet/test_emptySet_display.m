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

res = false;

% 0-dimensional empty set
O = emptySet(0)

% n-dimensional empty set
n = 2;
O = emptySet(2);

% only has to run through without errors...
res = true;

% ------------------------------ END OF CODE ------------------------------
