function S_out = plus(P,S)
% plus - overloaded '+' operator for the Minkowski addition of a polytope
%    and another set or vector
%
% Syntax:
%    S_out = P + S
%    S_out = plus(P,S)
%
% Inputs:
%    P - polytope object, numeric
%    S - contSet object, numeric
%
% Outputs:
%    S_out - polytope after Minkowski addition
%
% Example:
%    A = [2 1; -1 1; -2 -3; 0 -4; 2 -1]; b = ones(5,1);
%    P = polytope(A,b);
%    v = [2;1];
%    res = P + v;
%
%    figure; hold on;
%    plot(P);
%    plot(res,[1,2],'r');
%
% Reference:
%    [1] M. Wetzlinger, V. Kotsev, A. Kulmburg, M. Althoff. "Implementation
%        of Polyhedral Operations in CORA 2024", ARCH'24.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Viktor Kotsev
% Written:       20-June-2022
% Last update:   25-October-2023
% Last revision: 14-July-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[P,S] = reorderNumeric(P,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < P.precedence
    S_out = S + P;
    return
end

% check dimensions
equalDimCheck(P,S);

% dimension
n = dim(P);

% set tolerance for set checks
tol = 1e-10;

% check for fullspace
if representsa_(P,'fullspace',tol) || representsa_(S,'fullspace',tol)
    S_out = polytope.Inf(n); return
end

% polytope + polytope
if isa(S,'polytope')

    % check for empty set
    if representsa_(P,'emptySet',tol) || representsa_(S,'emptySet',tol)
        S_out = polytope.empty(n); return
    end
    % check if one of the sets is the origin
    if representsa_(P,'origin',tol)
        S_out = S; return
    elseif representsa_(S,'origin',tol)
        S_out = P; return
    end

    % check which representations are given
    if P.isHRep.val && S.isHRep.val
        S_out = aux_plus_Hpoly_Hpoly(P,S,n);
    elseif P.isVRep.val && S.isVRep.val
        S_out = aux_plus_Vpoly_Vpoly(P,S,n);
    else
        % convert to H-rep (since probably more useful for subsequent
        % operations) and compute sum
        constraints(P);
        constraints(S);
        S_out = aux_plus_Hpoly_Hpoly(P,S,n);
    end
    S_out = aux_setproperties(S_out,P,S);
    return
end

if isnumeric(S) && iscolumn(S)
    S_out = aux_plus_poly_point(P,S);
    S_out = aux_setproperties(S_out,P,S);
    return
end

% other set representations: convert to polytope
if isa(S,'zonotope') || isa(S,'interval') || isa(S,'conZonotope') || isa(S,'zonoBundle')
    S = polytope(S);
    S_out = P + S;
    return
end

throw(CORAerror('CORA:noops',P,S));

end


% Auxiliary functions -----------------------------------------------------

function P_out = aux_plus_Hpoly_Hpoly(P1,P2,n)

% all-zero matrices
PZ = zeros(size(P1.A_.val,1),size(P2.A_.val,2));
PZe = zeros(size(P1.Ae_.val,1),size(P2.Ae_.val,2));

% lift inequalities and equalities to higher-dimensional space
A = [P2.A_.val -P2.A_.val; PZ P1.A_.val];
b = [P2.b_.val; P1.b_.val];
Ae = [P2.Ae_.val -P2.Ae_.val; PZe P1.Ae_.val];
be = [P2.be_.val; P1.be_.val];

% project resulting polytope onto original dimensions
P_highdim = polytope(A,b,Ae,be);
P_out = project(P_highdim,1:n);

end

function P_out = aux_plus_Vpoly_Vpoly(P1,P2,n)
% Minkowski sum of two V-polytopes according to [1, (25)]

% read out vertices
V1 = P1.V_.val;
V2 = P2.V_.val;
numV1 = size(V1,2); numV2 = size(V2,2);

% init resulting vertices
V = zeros(n,numV1*numV2);

% add each vertex to each vertex
for i = 1:numV2
	V(:,((i-1)*numV1+1):i*numV1,:) = bsxfun(@plus, V1, V2(:,i));
end

% init resulting polytope
P_out = polytope(V);

end

function P_out = aux_plus_poly_point(P,S)
% Minkowski sum of a polytope and a vector

% avoid addition with matrices
if size(S,2) > 1
    throw(CORAerror('CORA:noops',P,S));
elseif size(S,1) ~= dim(P)
    throw(CORAerror('CORA:notSupported',...
        'Minkowski addition with scalar is not supported unless the set is 1-dimensional.'));
end

% copy polytope (and properties!)
P_out = polytope(P);

% shift offsets
if P_out.isHRep.val
    [~,P_out.b_.val,~,P_out.be_.val] = priv_plus_minus_vector(P.A_.val,P.b_.val,P.Ae_.val,P.be_.val,S);
end

% shift vertices if V representation is given
if P_out.isVRep.val
    P_out.V_.val = P.V_.val + S;
end

end

function P_out = aux_setproperties(P_out,P,S)
% infer set properties following [1, Table 1]

% currently only supported if S is a polytope
if isa(S,'polytope')

    % if both polytopes are bounded, then sum is also bounded
    if (~isempty(P.bounded.val) && P.bounded.val) ...
        && (~isempty(S.bounded.val) && S.bounded.val)
        P_out.bounded.val = true;
    end

    % If one of the polytopes is unbounded, then sum is also unbounded
    if (~isempty(P.bounded.val) && ~P.bounded.val) ...
        || (~isempty(S.bounded.val) && ~S.bounded.val)
        P_out.bounded.val = false;
    end

    % If one of the polytopes is non-degenerate, the sum is, too
    if (~isempty(P.fullDim.val) && P.fullDim.val) ...
        || (~isempty(S.fullDim.val) && S.fullDim.val)
        P_out.fullDim.val = true;
    end

end

end

% ------------------------------ END OF CODE ------------------------------
