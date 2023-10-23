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

% Authors:       Mark Wetzlinger
% Written:       30-November-2022
% Last update:   ---
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


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
