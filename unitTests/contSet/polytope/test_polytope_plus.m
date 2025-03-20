function res = test_polytope_plus
% test_polytope_plus - unit test function of Minkowski addition
%
% Syntax:
%    res = test_polytope_plus
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

% Authors:       Mark Wetzlinger, Viktor Kotsev
% Written:       30-November-2022
% Last update:   25-October-2023
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-14;

% 1D, fully empty + vector
P = polytope.empty(1);
z = 1;
P_sum = P + z;
assert(representsa(P_sum,'emptySet'));

% 1D, bounded + scalar
A = [1; -1]; b = [4; -2];
P = polytope(A,b);
z = 1;
P_sum = P + z;
P_true = polytope([1;-1],[5;-3]);
assert(isequal(P_sum,P_true,tol));


% 2D, bounded + vector (vertex representation)
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);
z = [2; 1];
P_sum = P + z;
V_sum = V + z;
P_true = polytope(V_sum);
assert(isequal(P_sum,P_true,tol));

% 2D, bounded + bounded (vertex representation)
V = [1 1; -1 1; -1 -1; 1 -1]';
P1 = polytope(V);
V = [0 1; 1 0; -1 0]';
P2 = polytope(V);
P_sum = compact(P1 + P2,'V');
V_true = [-1 2; 1 2; 2 1; 2 -1; -2 -1; -2 1]';
P_true = polytope(V_true);
assert(compareMatrices(P_sum.V,P_true.V));

% 2D, unbounded + bounded
A = [1 0; -1 0; 0 1]; b = [1;1;1];
P1 = polytope(A,b);
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;1];
P2 = polytope(A,b);
P_sum = P1 + P2;
A = [1 0; -1 0; 0 1]; b = [2;2;2];
P_true = polytope(A,b);
assert(P_true == P_sum);

% 2D, degenerate + degenerate
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;-1];
P1 = polytope(A,b);
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;-1];
P2 = polytope(A,b);
P_sum = P1 + P2;
A_true = [1 0; -1 0; 0 1; 0 -1]; b_true = [2;2;2;-2];
P_true = polytope(A_true,b_true);
assert(isequal(P_sum,P_true,tol));

% 2D, bounded + degenerate
A = [1 0; -1 0; 0 1; 0 -1]; b = 0.1*ones(4,1);
P1 = polytope(A,b);
A = [-1 0; 1 0]; b = [0.05;0.05]; Ae = [1 -1]; be = 0;
P2 = polytope(A,b,Ae,be);
P_sum = P1 + P2;
V_true = [-0.15 -0.15; -0.15 0.05; 0.05 -0.15; 0.05 0.05; ...
          -0.05 -0.05; -0.05 0.15; 0.15 -0.05; 0.15 0.15]';
P_true = polytope(V_true);
assert(isequal(P_sum,P_true,tol));


% 2D, bounded + interval
A = [2 1; -1 1; -2 -3; 0 -4; 2 -1]; b = ones(5,1);
P = polytope(A,b);
lb = [-2;1]; ub = [3;2];
I = interval(lb,ub);
P_sum = P + I;
assert(isa(P_sum,"polytope"));

% 2D, polytope + zonotope
A = [2 1; -1 1; -2 -3; 0 -4; 2 -1]; b = ones(5,1);
P = polytope(A,b);
isBounded(P);
c = [1;0]; G = [1 0; 0 1];
Z = zonotope(c,G);
P_sum = P + Z;
assert(isa(P_sum,"polytope") && (~isempty(P_sum.bounded.val) && P_sum.bounded.val));

% 2D, polytope + polyZonotope
A = [2 1; -1 1; -2 -3; 0 -4; 2 -1]; b = ones(5,1);
P = polytope(A,b);
c = [1;0]; G = [1 0; 0 1];
pZ = polyZonotope(c,G);
P_sum = P + pZ;
assert(isa(P_sum,"polyZonotope"));

% 2D, polytope + point
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
p = [1;2];
P_sum = P + p;
V_true = [2 4; 2 0; 0 2]';
P_true = polytope(V_true);
assert(isequal(P_sum,P_true,tol));


% 2D, addition with scalar
A = [0.3536 0; 0.5 -1; -0.3536 0; -0.5 1];
b = [sqrt(2); 1; 0; 1];
P = polytope(A,b);
assertThrowsAs(@plus,'CORA:notSupported',P,0);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
