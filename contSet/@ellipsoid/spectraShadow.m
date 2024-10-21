function SpS = spectraShadow(E)
% spectraShadow - converts an ellipsoid to a spectrahedral shadow
%
% Syntax:
%    SpS = spectraShadow(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    SpS - spectraShadow object
%
% Example:
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    SpS = spectraShadow(E);
%
% References:
%    [1] Ramana and Goldman, "Some Geometric Results in Semidefinite
%        Programming", Journal of Global Optimization, 1995.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       01-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% We need to 'separate' the inverse of the shape matrix first. How this can
% be done depends on the degeneracy of the ellipsoid:
if isFullDim(E)
    [A0, Ai] = aux_nondegenerate(E);
else
    [A0, Ai] = aux_degenerate(E);
end

% concatenate everything
A = [A0, cat(2,Ai{:})];

% instantiate spectraShadow
SpS = spectraShadow(A);

% additional properties
SpS.bounded.val = true;
SpS.emptySet.val = representsa_(E,'emptySet',1e-10);
SpS.fullDim.val = isFullDim(E);
SpS.center.val = center(E);

end


% Auxiliary functions -----------------------------------------------------

function [A0, Ai] = aux_nondegenerate(E)

n = dim(E);
L = chol(inv(E.Q));
    
% According to [1, Section 1.4], the ellipsoid defined as
%     (x-q)' L'L (x-q) < 1
% can equivalently be defined through the linear matrix inequality
%     [1-q'L'Lq+2q'L'Lx   x'L']
%     [Lx                 I   ] > 0
% We thus have
%      [1-q'L'Lq  0']       [(2q'L'L)(i)  L(:,i)']
% A0 = [0       I ], Ai = [L(:,i)    0      ]
A0 = [1-E.q'*(L'*L)*E.q sparse(1,n); sparse(n,1) speye(n)];
Ai = cell(1,n);

temp = 2*E.q'*(L'*L);
for i = 1:n
    Ai{i} = [temp(i), L(:,i)'; L(:,i), sparse(n,n)];
end

end

function [A0, Ai] = aux_degenerate(E)

n = dim(E);
% First, let us make our life a bit easier by "rotating" Q such that it
% becomes a diagonal matrix:
[U,D,V] = svd(E.Q);
% In theory, U=V. But more importantly, we have
% U'E(Q,q) = E(U'QU,U'q) = E(D,U'q)

% Now, D might not be invertible, but that's no biggy. We now know in
% which axis-directions the ellipsoid E(D,U'q) is degenerate, namely
% those axis-directions in which the diagonal matrix D has 0-entries.
% So, suppose D = diag(D_r,0,...,0), with D_r invertible. As it
% happens, the Moore-Penrose pseudo-inverse D+ of D coincides with
% diag(D_r^{-1},0,...,0), so we have
% E(D,U'q) = {x+U'q | x'*D+*x < 1, x_{r+1}=...=x_{n}=0}
% where r is the rank of D.
% This finally gives us the formula
% E(Q,q) = U*E(D,U'q) = {U*x+q | x'*D+*x < 1, x_{r+1}=...=x_{n}=0}
%        = {U*x+q | x'*D+*x < 1, [0 I_{n-r}]x=0}
%        = {y | (y-q)UD+U'(y-q) < 1, [0 I_{n-r}]U'y=[0 I_{n-r}]q}
%        = {y | (y-q)Q+(y-q) < 1, [0 I_{n-r}]U'y=[0 I_{n-r}]q}
% where the last equality follows from the very useful coincidence that
% Q+ = UD+U' (always the case for the svd).
% Therefore, we get something that is very close to the non-degenerate
% case, with a few equality conditions on top.

% We first deal with the 'usual' elliptic constraints:
% Instead of the Cholesky decomposition, we need to cheat a bit:
[U_pinv, D_pinv, V_pinv] = svd(pinv(E.Q));
% Again, in theory U_pinv = V_pinv
% So pinv(E.Q) = V_pinv * D_pinv * V_pinv'
% = V_pinv * sqrt(D_pinv) * sqrt(D_pinv) * V_pinv'
% So we can choose L = sqrt(D_pinv) * V_pinv'
L = sqrt(D_pinv) * V_pinv';

A0 = [1-E.q'*(L'*L)*E.q, sparse(1,n); sparse(n,1), speye(n)];
Ai = cell([1 n]);

temp = 2*E.q'*(L'*L);
for i = 1:n
    Ai{i} = [temp(i) L(:,i)'; L(:,i) sparse(n,n)];
end

% We now add the equality constraints:
r = rank(E.Q);
T = blkdiag(sparse(r,r),speye(n-r,n-r)) * U';
t = blkdiag(sparse(r,r),speye(n-r,n-r)) * E.q;

for i = 1:n
    A0 = blkdiag(A0,sparse([t(i) 0;0 -t(i)]));
    for j=1:n
        Ai{i} = blkdiag(Ai{i}, sparse([-T(i,j) 0;0 T(i,j)]));
    end
end

end

% ------------------------------ END OF CODE ------------------------------
