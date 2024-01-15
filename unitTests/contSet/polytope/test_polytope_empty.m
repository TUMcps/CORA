function res = test_polytope_empty
% test_polytope_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_polytope_empty
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
% Written:       16-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D
n = 1;
P = polytope.empty(n);
res(end+1,1) = representsa(P,'emptySet') && dim(P) == 1;

% 5D
n = 5;
P = polytope.empty(n);
res(end+1,1) = representsa(P,'emptySet') && dim(P) == 5;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
