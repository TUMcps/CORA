function [c_,G_,A_,b_,empty] = priv_conZonotope_supportFunc(A,b,Ae,be,B)
% priv_conZonotope_supportFunc - converts a polytope to a constrained
%    zonotope using support function evaluations
%
% Syntax:
%    [A,b,Ae,be,empty] = priv_conZonotope_supportFunc(A,b,Ae,be,B)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    B - (any) enclosure of polytope
%
% Outputs:
%    c_ - center of conZonotope
%    G_ - generator matrix of conZonotope
%    A_ - constraint matrix of conZonotope
%    b_ - constraint offset of conZonotope
%    empty - true/false whether polytope is the empty set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute center and generators
c_ = center(B);
G = diag(0.5 * (supremum(B) - infimum(B)));

% number of constraints
nrIneq = numel(b);
nrEq = numel(be);

% read out all constraints
A_all = [A; Ae];
b_all = [b; be];

% compute lower bound in the direction of halfspaces
sigma = zeros(nrIneq,1);
for a=1:nrIneq
    sigma(a) = supportFunc_(B,A_all(a,:)','lower');
    % any lower bound is Inf -> polytope is empty
    if sigma(a) == Inf
        c_ = []; G_ = []; A_ = []; b_ = [];
        empty = true;
        return
    end
end
% same for equality constraints (no need to compute the value)
sigma = [sigma; be];

% construct constrained zonotope object according to eq. (21) in [1]
G_ = [G, zeros(size(G,1),nrIneq+nrEq)];
A_ = [A_all*G, diag((sigma-b_all)./2)];
b_ = (b_all+sigma)./2 - A_all*c_;

% non-empty
empty = false;

% ------------------------------ END OF CODE ------------------------------
