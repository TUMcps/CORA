function res = test_zonotope_contains
% test_zonotope_contains - unit test function of contains
%
% Syntax:
%    res = test_zonotope_contains
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
% Written:       12-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty zonotope
% Z = zonotope.empty(2);
% res(end+1,1) = contains(Z,Z);

% 2D zonotopes
c = [0; 1]; G = [1 2; -1 0];
Z1 = zonotope(c,G);
c = [-1; 1.5]; G = [0.2 0; -0.1 0.1];
Z2 = zonotope(c,G);
res(end+1,1) = contains(Z1,Z2);
res(end+1,1) = ~contains(Z2,Z1);

% inner zonotope is just a point
Z2 = zonotope(c);
res(end+1,1) = contains(Z1,Z2);

% both zonotope are just points
Z1 = zonotope(zeros(4,1));
Z2 = zonotope(ones(4,1));
res(end+1,1) = contains(Z1,Z1);
res(end+1,1) = ~contains(Z1,Z2);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
