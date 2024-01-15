function res = test_zonotope_zonoBundle
% test_zonotope_zonoBundle - unit test function of conversion to zonoBundle
%
% Syntax:
%    res = test_zonotope_zonoBundle
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set
Z = zonotope.empty(2);
zB = zonoBundle(Z);
res = representsa(zB,'emptySet') && dim(zB) == 2;

% instantiate zonotope
c = [1;1;-1];
G = [2 -1 3 4; 0 2 3 -1; -1 0 0 2];
Z = zonotope(c,G);

% convert to zonoBundle
zB = zonoBundle(Z);

% check converted set
res(end+1,1) = zB.parallelSets == 1;
res(end+1,1) = isequal(zB.Z{1},Z);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
