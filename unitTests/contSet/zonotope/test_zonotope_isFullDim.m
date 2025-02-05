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

% Authors:       Mark Wetzlinger, Adrian Kulmburg
% Written:       27-July-2021
% Last update:   04-February-2024 (AK, added degenerate cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% empty case
Z = zonotope.empty(2);
[res_subspace, subspace] = isFullDim(Z);
assert(~isFullDim(Z));
assert(~res_subspace);
assert(isempty(subspace));

% full-rank matrix
c = [0; 0];
G = [1 2; 2 1];
Z = zonotope(c,G);
[res_subspace, subspace] = isFullDim(Z);
assert(isFullDim(Z));
assert(res_subspace);
assert(size(subspace,2)==2);

% degenerate zonotope
G = [1 0; 0 0];
Z = zonotope(c,G);
[res_subspace, subspace] = isFullDim(Z);
true_subspace = [1;0];
same_subspace = rank([subspace true_subspace], 1e-6) == size(subspace,2);
assert(~isFullDim(Z));
assert(~res_subspace);
assert(same_subspace);

% almost degenerate zonotope
eps = 1e-8;
c = [0;0];
G = [1 1-eps; 1 1];
Z = zonotope(c,G);
assert(~isFullDim(Z));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
