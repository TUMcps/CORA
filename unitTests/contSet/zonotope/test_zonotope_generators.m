function res = test_zonotope_generators
% test_zonotope_generators - unit test function of generators
%
% Syntax:
%    res = test_zonotope_generators
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

% empty zonotope
Z = zonotope.empty(2);
G = generators(Z);
assert(isempty(G) && isnumeric(G) && all(size(G) == [2,0]));

% 2D zonotope
c = [-2; 1];
G = [2 4 5 3 3; 0 3 5 2 3];
Z = zonotope(c,G);
G_ = generators(Z);
assert(compareMatrices(G,G_));

% zonotope without generators
c = [1;0;0;-2];
Z = zonotope(c);
G_ = generators(Z);
assert(all(size(G_) == [4,0]));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
