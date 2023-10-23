function res = test_fullspace_display
% test_fullspace_display - unit test function of display
%
% Syntax:
%    res = test_fullspace_display
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

% n-dimensional empty set
n = 2;
fs = fullspace(n)

% only has to run through without errors...
res = true;

% ------------------------------ END OF CODE ------------------------------
