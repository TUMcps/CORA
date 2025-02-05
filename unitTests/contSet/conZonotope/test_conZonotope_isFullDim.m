function res = test_conZonotope_isFullDim
% test_conZonotope_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = test_conZonotope_isFullDim
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
% Written:       21-April-2023
% Last update:   04-February-2024 (AK, added degenerate case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty conZonotope object
cZ = conZonotope.empty(2);
[res_subspace, subspace] = isFullDim(cZ);
assert(~isFullDim(cZ));
assert(~res_subspace);
assert(isempty(subspace))

% constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1]; b = 1;
cZ = conZonotope(Z,A,b);
[res_subspace, subspace] = isFullDim(cZ);
assert(isFullDim(cZ));
assert(res_subspace);
assert(size(subspace, 2) == 2);

% degenerate constrained zonotope
Z = [0 3 0 1; 1 0 0 0];
A = [1 0 1]; b = 1;
cZ = conZonotope(Z,A,b);
[res_subspace, subspace] = isFullDim(cZ);
true_subspace = [1;0];
same_subspace = rank([subspace true_subspace], 1e-6) == size(subspace,2);
assert(~isFullDim(cZ));
assert(~res_subspace);
assert(same_subspace);

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
