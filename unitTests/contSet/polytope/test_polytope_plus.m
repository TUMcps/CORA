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

res = [];

% init polytope
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);

% translate
z = [2; 1];
V_ = V + z;
P_ = polytope(V_);

% plot polytopes
% figure; hold on;
% plot(P); plot(P_,[1,2],'r');

% constraint matrices should be equal
res(end+1,1) = compareMatrices(P.A' ./ vecnorm(P.A'),P_.A' ./ vecnorm(P_.A'),1e-14);

% comparison to manual translation
val1 = zeros(size(P.A))';
val2 = zeros(size(P.A))';
for i=1:size(P.A,1)
    val1(:,i) = (P.A(i,:) * P_.b(i))';
    val2(:,i) = (P_.A(i,:) * (P.b(i) + P.A(i,:)*z))';
end
res(end+1,1) = compareMatrices(val1,val2,1e-14);

% two vertex representations
P1 = polytope([1 1; -1 1; -1 -1; 1 -1]');
P2 = polytope([0 1; 1 0; -1 0]');
P_sum = compact(P1 + P2,'V');
P_true = polytope([-1 2; 1 2; 2 1; 2 -1; -2 -1; -2 1]');
res(end+1,1) = compareMatrices(P_sum.V.val,P_true.V.val);

% unbounded case
P1 = polytope([1 0; -1 0; 0 1], [1;1;1]);
P2 = polytope([1 0; -1 0; 0 1;0 -1], [1;1;1;1]);

P_sum = P1 + P2;
P_true = polytope([1 0; -1 0; 0 1], [2;2;2]);
res(end+1,1) = P_true == P_sum;

% degenerate case
P1 = polytope([1 0; -1 0; 0 1; 0 -1], [1;1;1;-1]);
P2 = polytope([1 0; -1 0; 0 1; 0 -1], [1;1;1;-1]);

P_sum = P1 + P2;
P_true = polytope([1 0; -1 0; 0 1; 0 -1], [2;2;2;-2]);
res(end+1,1) = P_sum == P_true;

% degenerate case
P1 = polytope([1 0; -1 0; 0 1; 0 -1],0.1*ones(4,1));
P2 = polytope([-1 0; 1 0],[0.05;0.05],[1 -1],0);

P_sum = P1 + P2;
P_true = polytope([-0.15 -0.15; -0.15 0.05; 0.05 -0.15; 0.05 0.05; ...
        -0.05 -0.05; -0.05 0.15; 0.15 -0.05; 0.15 0.15]');
res(end+1,1) = P_sum == P_true;


% polytope + interval
P = polytope([2 1; -1 1; -2 -3; 0 -4; 2 -1], ones(5,1));
I =  interval([-2;1],[3;2]);
PI = P + I;
res(end+1,1) = isa(PI,"polytope");

% polytope + zonotope
P = polytope([2 1; -1 1; -2 -3; 0 -4; 2 -1], ones(5,1));
isBounded(P);
Z = zonotope([1 1 0; 0 0 1]);
PZ = P + Z;
res(end+1,1) = isa(PZ,"polytope") && (~isempty(PZ.bounded.val) && PZ.bounded.val);

% polytope + polyZonotope
P = polytope([2 1; -1 1; -2 -3; 0 -4; 2 -1], ones(5,1));
Z = polyZonotope([1;0], [1 0; 0 1]);
PZ = P + Z;
res(end+1,1) = isa(PZ,"polyZonotope");

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
