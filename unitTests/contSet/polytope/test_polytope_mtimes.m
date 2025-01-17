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

tol = 1e-14;

% 1D, vertex instantiation
V = [1,3];
P = polytope(V);
M = -2;
P_mtimes = M*P;
V = vertices(P_mtimes);
V_true = [-2,-6];
assert(compareMatrices(V,V_true,tol));


% 2D, vertex instantiation
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);
M = [2 1; -1 0];
P_mtimes = M*P;
V_mapped = vertices(P_mtimes);
V_true = M*V;
% compare vertices
assert(compareMatrices(V_true,V_mapped,tol));

% 2D, vertex instantiation, projection
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);
M = [2 1];
P_mtimes = M*P;
V_mapped = vertices(P_mtimes);
V_true = [min(M*V), max(M*V)]; % because 1D
% compare vertices
assert(compareMatrices(V_true,V_mapped,tol));

% 2D, vertex instantiation, lifting
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);
M = [2 1; -1 0; 3 -1];
P_mtimes = M*P;
V_mapped = vertices(P_mtimes);
V_true = M*V;
% compare vertices
assert(compareMatrices(V_true,V_mapped,tol));


% 2D, unbounded
A = [1 0; -1 0; 0 1]; b = [1;1;1];
P = polytope(A,b);
M = [2 1; -1 0];
P_mtimes = M * P;
assert(~isBounded(P_mtimes));

% 2D, degenerate
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;-1];
P = polytope(A,b);
M = [2 1; -1 0];
P_mtimes = M * P;
assert(~isFullDim(P_mtimes));

% 2D, box, scaling matrix
A = [1 0;-1 0; 0 1; 0 -1]; b = [2;2;2;2];
P = polytope(A,b);
M = [2 0; 0 4];
P_mtimes = mtimes(M,P);
A_true = [1 0;-1 0; 0 1; 0 -1]; b_true = [4;4;8;8];
P_true = polytope(A_true,b_true);
assert(P_mtimes == P_true);

% check multiplication with scalar (left and right)
P = polytope([1 0; -1 1; -1 -1], [1; 1; 1]);
P_mtimes = 2*eye(2) * P;
P_mtimes_scalar = 2 * P;
assert(P_mtimes == P_mtimes_scalar);

P_mtimes_scalar = P * 2;
assert(P_mtimes == P_mtimes_scalar);

% 2D, lifting to 3D
Z = zonotope(zeros(2,1),[1 1 0; 1 2 1]);
P = polytope(Z);
M = [1 -1; 0 1; 2 1];
P_mtimes = M * P;
P_true = polytope(M*Z);
assert(P_true == P_mtimes);

% 2D, degenerate, lifting to 3D
Z = zonotope(zeros(2,1),[1; -1]);
P = polytope(Z);
M = [1 -1; 0 1; 2 1];
P_mtimes = M * P;
P_true = polytope(M*Z);
assert(P_true == P_mtimes);


% 3D, map with all-zero matrix
P = polytope([1 0 0; 0 1 0; 0 0 1; -1 -1 -1],ones(4,1));
M = zeros(3);
P_mtimes = M * P;
P_true = polytope(zeros(3,1));
assert(P_true == P_mtimes);

% 3D, projection onto 2D
Z = zonotope(zeros(3,1),[1 1 0; 1 2 1; 0 1 1]);
P = polytope(Z);
M = [1 -1 0; 0 1 1];
P_mtimes = M * P;
P_true = polytope(M*Z);
assert(P_true == P_mtimes);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
