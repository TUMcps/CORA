function P_out = box(P)
% box - computes an enclosing axis-aligned box represented as a polytope in
%    halfspace representation
%
% Syntax:
%    P_out = box(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    P_out - polytope object 
%
% Example:
%    A = [1 2; -2 1; -2 -2; 3 -1];
%    b = ones(4,1);
%    P = polytope(A,b);
%    B = box(P);
%
%    figure; hold on;
%    plot(P);
%    plot(B,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       16-May-2022
% Last update:   14-December-2022 (MW, unbounded case, MOSEK support, remove equality constraints)
% Last revision: 12-July-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% fullspace case
if representsa_(P,'fullspace',0)
    P_out = polytope.Inf(dim(P)); return
end

% the computation of the box outer approximation is much faster for the
% vertex representation, so we first check for that
if P.isVRep.val
    % --- V representation: take the min and max along each dimension
    P_out = aux_box_V(P);
else
    % --- H representation: compute the support function in the direction
    % of all 2n plus/minus axis-aligned basis vectors
    P_out = aux_box_H(P);
end

end


% Auxiliary functions -----------------------------------------------------

function P_out = aux_box_V(P)
% computation of box enclosure for a polytope in vertex representation

% read out dimension
n = dim(P);

% compute lower and upper bound
ub = max(P.V,[],2);
lb = min(P.V,[],2);

% construct constraint matrix and offset
A = [eye(n); -eye(n)];
b = [ub; -lb];
P_out = polytope(A, b);

% set properties
P_out.emptySet.val = false; % cannot be empty, otherwise V would be empty
P_out.bounded.val = any(isinf(ub)) || any(isinf(lb));
% 1D polytopes support -Inf/Inf vertices...
if n == 1
    P_out.fullDim.val = max(P.V_.val) - min(P.V_.val) > 0;
else
    P_out.fullDim.val = rank(P.V_.val,1e-10) == n;
end

end

function P_out = aux_box_H(P)
% computation of box enclosure for a polytope in halfspace representation

% read out dimension
n = dim(P);

% init bounds
ub = Inf(n,1);
lb = -Inf(n,1);

% loop over all 2n positive/negative basis vectors
for i = 1:n
    % i-th basis vector
    e_i = [zeros(i-1,1);1;zeros(n-i,1)];
    % maximize
    ub(i) = supportFunc_(P,e_i,'upper');
    if ub(i) == -Inf
        % only enters here if the polytope is the empty set
        P.emptySet.val = true;
        P.bounded.val = true;
        P.fullDim.val = false;
        P.V_.val = zeros(n,0);
        P.isVRep.val = true;
        P.minVRep.val = true;
        % init return set
        P_out = polytope.empty(n);
        return
    end
    % minimize
    lb(i) = supportFunc_(P,e_i,'lower');
end

% construct output arguments
A = [eye(n); -eye(n)];
b = [ub; -lb];

% remove unbounded constraints
bounded = ~isinf(b);

% rewrite properties
A = A(bounded,:);
b = b(bounded);

% remove equality constraints
Ae = [];
be = [];

% set properties: input polytope
P.emptySet.val = false;
P.bounded.val = all(bounded);

% construct box
P_out = polytope(A,b,Ae,be);
% set properties: output polytope
P_out.minHRep.val = true;
P_out.emptySet.val = false;
P_out.bounded.val = all(bounded);
P_out.fullDim.val = ~any(withinTol(lb,ub,1e-10));

end

% ------------------------------ END OF CODE ------------------------------
