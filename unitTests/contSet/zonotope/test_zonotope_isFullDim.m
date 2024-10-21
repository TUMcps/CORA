function res = test_zonotope_isFullDim
% test_zonotope_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = test_zonotope_isFullDim
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
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% empty case
Z = zonotope.empty(2);
assert(~isFullDim(Z));

% full-rank matrix
c = [0; 0];
G = [1 2; 2 1];
Z = zonotope(c,G);
assert(isFullDim(Z));

% degenerate zonotope
G = [1 0; 0 0];
Z = zonotope(c,G);
assert(~isFullDim(Z));

% almost degenerate zonotope
eps = 1e-8;
c = [0;0];
G = [1 1-eps; 1 1];
Z = zonotope(c,G);
assert(~isFullDim(Z));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
