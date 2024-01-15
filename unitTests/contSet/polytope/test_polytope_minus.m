function res = test_polytope_minus
% test_polytope_minus - unit test function of minus
%
% Syntax:
%    res = test_polytope_minus
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

% 2D, vertex instantiation
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);
% translate by vector
z = [2; 1];
V_ = V - z;
P_ = polytope(V_);

% 1. check: constraint matrices should be equal
res(end+1,1) = compareMatrices(P.A' ./ vecnorm(P.A'),...
    P_.A' ./ vecnorm(P_.A'),1e-14);

% 2. check: comparison to manual translation
val1 = zeros(size(P.A))';
val2 = zeros(size(P.A))';
for i=1:size(P.A,1)
    val1(:,i) = (P.A(i,:) * P_.b(i))';
    val2(:,i) = (P_.A(i,:) * (P.b(i) - P.A(i,:)*z))';
end
res(end+1,1) = compareMatrices(val1,val2,1e-14);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
