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

res = true(0);
tol = 1e-14;

% 1D, fully empty + vector
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
z = 1;
P_sum = P + z;
res(end+1,1) = isemptyobject(P_sum);


% 2D, bounded + vector (vertex representation)
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);
z = [2; 1];
V_sum = V + z;
P_sum = polytope(V_sum);
% constraint matrices should be equal
res(end+1,1) = compareMatrices(P.A' ./ vecnorm(P.A'),...
    P_sum.A' ./ vecnorm(P_sum.A'),tol);
% comparison to manual translation
val1 = zeros(size(P.A))';
val2 = zeros(size(P.A))';
for i=1:size(P.A,1)
    val1(:,i) = (P.A(i,:) * P_sum.b(i))';
    val2(:,i) = (P_sum.A(i,:) * (P.b(i) + P.A(i,:)*z))';
end
res(end+1,1) = compareMatrices(val1,val2,tol);

% 2D, bounded + bounded (vertex representation)
V = [1 1; -1 1; -1 -1; 1 -1]';
P1 = polytope(V);
V = [0 1; 1 0; -1 0]';
P2 = polytope(V);
P_sum = compact(P1 + P2,'V');
V_true = [-1 2; 1 2; 2 1; 2 -1; -2 -1; -2 1]';
P_true = polytope(V_true);
res(end+1,1) = compareMatrices(P_sum.V.val,P_true.V.val);

% 2D, unbounded + bounded
A = [1 0; -1 0; 0 1]; b = [1;1;1];
P1 = polytope(A,b);
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;1];
P2 = polytope(A,b);
P_sum = P1 + P2;
A = [1 0; -1 0; 0 1]; b = [2;2;2];
P_true = polytope(A,b);
res(end+1,1) = P_true == P_sum;

% 2D, degenerate + degenerate
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;-1];
P1 = polytope(A,b);
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;-1];
P2 = polytope(A,b);
P_sum = P1 + P2;
A_true = [1 0; -1 0; 0 1; 0 -1]; b_true = [2;2;2;-2];
P_true = polytope(A_true,b_true);
res(end+1,1) = P_sum == P_true;

% 2D, bounded + degenerate
A = [1 0; -1 0; 0 1; 0 -1]; b = 0.1*ones(4,1);
P1 = polytope(A,b);
A = [-1 0; 1 0]; b = [0.05;0.05]; Ae = [1 -1]; be = 0;
P2 = polytope(A,b,Ae,be);
P_sum = P1 + P2;
V_true = [-0.15 -0.15; -0.15 0.05; 0.05 -0.15; 0.05 0.05; ...
          -0.05 -0.05; -0.05 0.15; 0.15 -0.05; 0.15 0.15]';
P_true = polytope(V_true);
res(end+1,1) = P_sum == P_true;


% 2D, bounded + interval
A = [2 1; -1 1; -2 -3; 0 -4; 2 -1]; b = ones(5,1);
P = polytope(A,b);
lb = [-2;1]; ub = [3;2];
I = interval(lb,ub);
P_sum = P + I;
res(end+1,1) = isa(P_sum,"polytope");

% 2D, polytope + zonotope
A = [2 1; -1 1; -2 -3; 0 -4; 2 -1]; b = ones(5,1);
P = polytope(A,b);
isBounded(P);
c = [1;0]; G = [1 0; 0 1];
Z = zonotope(c,G);
P_sum = P + Z;
res(end+1,1) = isa(P_sum,"polytope") ...
    && (~isempty(P_sum.bounded.val) && P_sum.bounded.val);

% 2D, polytope + polyZonotope
A = [2 1; -1 1; -2 -3; 0 -4; 2 -1]; b = ones(5,1);
P = polytope(A,b);
c = [1;0]; G = [1 0; 0 1];
pZ = polyZonotope(c,G);
P_sum = P + pZ;
res(end+1,1) = isa(P_sum,"polyZonotope");


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
