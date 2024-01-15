function res = test_polytope_mtimes
% test_polytope_mtimes - unit test function of linear map
%
% Syntax:
%    res = test_polytope_mtimes
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
% Written:       30-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-14;

% 1D, vertex instantiation
V = [1,3];
P = polytope(V);
M = -2;
P_ = M*P;
V = vertices(P_);
V_true = [-2,-6];
res(end+1,1) = compareMatrices(V,V_true,tol);


% 2D, vertex instantiation
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);
M = [2 1; -1 0];
P_mapped = M*P;
V_mapped = vertices(P_mapped);
V_ = M*V;
P_ = polytope(V_);
% compare vertices
res(end+1,1) = compareMatrices(V_,V_mapped,tol);
% normalize constraints and compare
Ab_mapped = [P_mapped.A, P_mapped.b]';
Ab_mapped = Ab_mapped ./ vecnorm(Ab_mapped);
Ab_ = [P_.A, P_.b]';
Ab_ = Ab_ ./ vecnorm(Ab_);
res(end+1,1) = compareMatrices(Ab_mapped,Ab_,tol);

% 2D, unbounded
A = [1 0; -1 0; 0 1]; b = [1;1;1];
P = polytope(A,b);
M = [2 1; -1 0];
P_ = M * P;
res(end+1,1) = ~isBounded(P_);

% 2D, degenerate
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;-1];
P = polytope(A,b);
M = [2 1; -1 0];
P_ = M * P;
res(end+1,1) = ~isFullDim(P_);

% 2D, box, scaling matrix
A = [1 0;-1 0; 0 1; 0 -1]; b = [2;2;2;2];
P = polytope(A,b);
M = [2 0; 0 4];
P_ = mtimes(M,P);
A_true = [1 0;-1 0; 0 1; 0 -1]; b_true = [4;4;8;8];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_ == P_true;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
