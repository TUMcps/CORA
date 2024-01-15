function res = test_polytope_minkDiff
% test_polytope_minkDiff - unit test function of minkDiff
%
% Syntax:
%    res = test_polytope_minkDiff
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
% Written:       07-September-2022
% Last update:   01-December-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-10;

% 1D, bounded - bounded
A = [1; -1]; b = [1; 1];
P1 = polytope(A,b);
A = [1; -1]; b = [0.2; 0.7];
P2 = polytope(A,b);
P_minkDiff = minkDiff(P1,P2);
A_true = [1;-1]; b_true = [0.8; 0.3];
P_true = polytope(A_true,b_true);
res(end+1,1) = isequal(P_minkDiff,P_true,tol);

% 1D, fullspace - unbounded
A = zeros(0,1); b = zeros(0,0);
P1 = polytope(A,b);
A = 1; b = 1;
P2 = polytope(A,b);
P_minkDiff = minkDiff(P1,P2);
res(end+1,1) = isequal(P_minkDiff,P1,tol);

% 1D, fullspace - fullspace
A = zeros(0,1); b = zeros(0,0);
P1 = polytope(A,b);
P_minkDiff = minkDiff(P1,P1);
res(end+1,1) = isequal(P_minkDiff,P1,tol);


% 2D, origin - origin
V = zeros(2,1);
P1 = polytope(V);
P_minkDiff = minkDiff(P1,P1);
res(end+1,1) = isequal(P_minkDiff,P1,tol);

% 2D, set - origin
A = [1 1; -2 1; 0 -2]; b = [1;1;1];
P1 = polytope(A,b);
V = zeros(2,1);
P2 = polytope(V);
P_minkDiff = minkDiff(P1,P2);
res(end+1,1) = isequal(P_minkDiff,P1,tol);

% 2D, convert zonotope to polytopes
c = [1; 1]; G = [1 0 1; 0 1 1];
Z1 = zonotope(c,G);
P1 = polytope(Z1);
c = [-0.5; 1]; G = [0.5 0; 0 0.5];
Z2 = zonotope(c,G);
P2 = polytope(Z2);
% compute Minkowski difference
P_minkDiff = minkDiff(P1,P2);
% true result in zonotope representation and converted to a polytope
c_true = [1.5; 0]; G_true = [0.5 0 1; 0 0.5 1];
Z_true = zonotope(c_true,G_true);
P_true = polytope(Z_true);
res(end+1,1) = P_minkDiff == P_true;

% 2D, bounded - unbounded
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;1];
P1 = polytope(A,b);
A = [1 0; -1 0; 0 1]; b = [1;1;1];
P2 = polytope(A,b);
P_minkDiff = minkDiff(P1,P2);
res(end+1,1) = representsa(P_minkDiff,'emptySet');

% 2D, degenerate - unbounded
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;0;0];
P1 = polytope(A,b);
A = [1 0; -1 0; 0 1]; b = [1;1;1];
P2 = polytope(A,b);
P_minkDiff = minkDiff(P1,P2);
res(end+1,1) = representsa(P_minkDiff,'emptySet');

% 2D, unbounded - degenerate
A = [1 0; -1 0; 0 1]; b = [1;1;1];
P1 = polytope(A,b);
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;0;0];
P2 = polytope(A,b);
P_minkDiff = minkDiff(P1,P2);
A_true = [1 0; -1 0; 0 1]; b_true = [0;0;1];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_minkDiff == P_true;

% 2D, unbounded - unbounded
A = [1 0; -1 0; 0 -1]; b = [1;1;1];
P1 = polytope(A,b);
A = [1 0; -1 0; 0 1]; b = [1;1;1];
P2 = polytope(A,b);
P_minkDiff = minkDiff(P1,P2);
res(end+1,1) = representsa(P_minkDiff, 'emptySet');

% 2D, two identical polytopes, unbounded
A = [1 0; -1 0; 0 1]; b = [1;1;1];
P1 = polytope(A,b);
P_minkDiff = minkDiff(P1,P1);
A = [1 0; -1 0; 0 1]; b = [0;0;0];
P_true = polytope(A,b);
res(end+1,1) = P_minkDiff == P_true;

% 2D, two identical polytopes, bounded
A = [-1 0; 0 -1]; b = [0;0];
P1 = polytope(A,b);
P_minkDiff = minkDiff(P1,P1);
res(end+1,1) = P_minkDiff == P1;

% 2D, degenerate - degenerate (equality constraints)
A = [1 0 0; -1 0 0; 0 1 0; 0 -1 0]; b = [1;1;1;1]; Ae = [0 0 1]; be = 0;
P1 = polytope(A,b,Ae,be);
A = [1 0 0; -1 0 0; 0 1 0; 0 -1 0]; b = 0.5*[1;1;1;1]; Ae = [0 0 1]; be = 0;
P2 = polytope(A,b,Ae,be);
P_minkDiff = minkDiff(P1,P2);
A_true = [1 0 0; -1 0 0; 0 1 0; 0 -1 0]; b_true = 0.5*[1;1;1;1];
Ae_true = [0 0 1]; be_true = 0;
P_true = polytope(A_true,b_true,Ae_true,be_true);
res(end+1,1) = P_minkDiff == P_true;

% 2D, degenerate - vector (only inequality constraints)
% ...line from (-1,0) to (1,0)
A = [1 0; -1 0]; b = [1;1]; Ae = [0 1]; be = 0;
P1 = polytope(A,b,Ae,be);
P2 = [-1;1];
P_minkDiff = minkDiff(P1,P2);
A_true = [1 0; -1 0]; b_true = [2;0]; Ae_true = [0 1]; be_true = -1;
P_true = polytope(A_true,b_true,Ae_true,be_true);
% ...line from (0,-1) to (2,-1)
res(end+1,1) = P_minkDiff == P_true;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
