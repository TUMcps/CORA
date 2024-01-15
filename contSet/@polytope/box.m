function P_out = box(P)
% box - computes an enclosing axis-aligned box
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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% fully empty object
if representsa_(P,'fullspace',0)
    P_out = polytope.Inf(dim(P)); return
end

% dimension
n = dim(P);
% check whether vertex representation is known
V = P.V.val;

% quick calculation using V-representation
if ~isempty(V)
    ub = max(V,[],2);
    lb = min(V,[],2);

    A = [eye(n); -eye(n)];
    b = [ub; -lb];
    P_out = polytope(A, b);

    % set properties
    P_out.emptySet.val = false; % cannot be empty, otherwise V would be empty
    P_out.bounded.val = any(isinf(ub)) && any(isinf(lb));
    P_out.fullDim.val = rank(V,1e-10) == n;
    return;
end

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
        % empty set
        P.emptySet.val = true;
        P.bounded.val = true;
        P.fullDim.val = false;
        P.V.val = zeros(n,0);
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

% ------------------------------ END OF CODE ------------------------------
