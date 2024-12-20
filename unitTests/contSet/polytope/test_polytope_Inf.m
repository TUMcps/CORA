function res = test_polytope_Inf
% test_polytope_Inf - unit test function of R^n instantiation
%
% Syntax:
%    res = test_polytope_Inf
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

% 1D
n = 1;
P = polytope.Inf(n);
assert(representsa(P,'fullspace') && dim(P) == 1);

% 5D
n = 5;
P = polytope.Inf(n);
assert(representsa(P,'fullspace') && dim(P) == 5);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
