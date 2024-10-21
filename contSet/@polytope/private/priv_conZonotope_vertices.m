function [c_,G_,A_,b_] = priv_conZonotope_vertices(A,b,Ae,be,V)
% priv_conZonotope_vertices - converts a polytope to a constrained
%    zonotope using vertex enumeration
%
% Syntax:
%    [A,b,Ae,be,empty] = priv_conZonotope_vertices(A,b,Ae,be,V)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    V - vertex representation
%
% Outputs:
%    c_ - center of conZonotope
%    G_ - generator matrix of conZonotope
%    A_ - constraint matrix of conZonotope
%    b_ - constraint offset of conZonotope
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

% number of constraints
nrIneq = numel(b);
nrEq = numel(be);

% read out all constraints
A_all = [A; Ae];
b_all = [b; be];

% calculate a bounding box for the constrained zonotope
minV = min(V,[],2);
maxV = max(V,[],2);

% compute center and generator matrix
c_ = 0.5 * (maxV + minV);
G = diag(0.5 * (maxV - minV));

% Calculate the lower bound sigma for A*x \in [sigma,b] (Thm. 1 in [1])
sigma = min(A_all*V,[],2);

% Construct constrained zonotope object according to eq. (21) in [1]
G_ = [G, zeros(size(G,1),nrIneq+nrEq)];
A_ = [A_all*G, diag((sigma-b_all)./2)];
b_ = (b_all+sigma)./2 - A_all*c_;

% ------------------------------ END OF CODE ------------------------------
