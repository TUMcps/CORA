function res = test_block_operation()
% test_block_operation - unit test function for pair-wise evaluation
%    of a set operation on two cell arrays of sets
%
% Syntax:
%    res = test_block_operation()
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
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-12;

% single pair of sets
S1 = interval([2;1],[3;1]);
S2 = interval([4;-2],[4;2]);

% evaluate different operations
S_res = block_operation(@plus,S1,S2);
S_true = S1 + S2;
assert(isequal(S_res,S_true,tol));
S_res = block_operation(@cartProd,S1,S2);
S_true = cartProd(S1,S2);
assert(isequal(S_res,S_true,tol));
S_res = block_operation(@convHull,S1,S2);
S_true = convHull(S1,S2);
assert(isequal(S_res,S_true,tol));

% multiple sets
S1 = {zonotope(1,0); zonotope([2;1],[1 0 -1; 1 1 2])};
S2 = {zonotope(-1,4); zonotope([-1;0],[-1 0; 1 2])};

% evaluate different operations
S_res = block_operation(@plus,S1,S2);
S_true = arrayfun(@(i) S1{i} + S2{i}, 1:2, 'UniformOutput', false);
assert(all(arrayfun(@(i) isequal(S_res{i},S_true{i},tol), 1:2, 'UniformOutput', true)));
S_res = block_operation(@cartProd,S1,S2);
S_true = arrayfun(@(i) cartProd(S1{i},S2{i}), 1:2, 'UniformOutput', false);
assert(all(arrayfun(@(i) isequal(S_res{i},S_true{i},tol), 1:2, 'UniformOutput', true)));
S_res = block_operation(@convHull,S1,S2);
S_true = arrayfun(@(i) convHull(S1{i},S2{i}), 1:2, 'UniformOutput', false);
assert(all(arrayfun(@(i) isequal(S_res{i},S_true{i},tol), 1:2, 'UniformOutput', true)));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
