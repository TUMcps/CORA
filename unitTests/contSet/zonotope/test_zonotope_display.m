function res = test_zonotope_display
% test_zonotope_display - unit test function of display
%
% Syntax:
%    res = test_zonotope_display
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
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% empty zonotope
n = 2;
Z = zonotope.empty(n)

% 2D zonotope
c = [-2; 1];
G = [2 4 5 3 3; 0 3 5 2 3];
Z = zonotope(c,G)

% no generator matrix
Z = zonotope(c)

% many generators
G = ones(2,25);
zonotope(c,G)

% ------------------------------ END OF CODE ------------------------------
