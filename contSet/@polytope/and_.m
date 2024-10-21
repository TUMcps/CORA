function P_out = and_(P,S,type,varargin)
% and_ - computes the intersection of a polytope and another set
%    note: the resulting representation is not necessarily minimal!
%
% Syntax:
%    P = and_(P,S,type)
%
% Inputs:
%    P - polytope object
%    S - contSet object or numerical vector
%    type - 'exact' or 'approx'
%
% Outputs:
%    P - polytope object
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([-1 0],-1,[1 1],2);
%
%    res = P1 & P2;
%
%    figure; hold on
%    xlim([-2,4]); ylim([-4,4]);
%    plot(P2,[1,2],'r','LineWidth',3);
%    plot(P1,[1,2],'b');
%    plot(res,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and, conZonotope/and_

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       09-May-2022
% Last update:   14-December-2022 (MW, bug fix, add equality constraints)
%                23-December-2023 (MW, support intersection with numeric)
% Last revision: 28-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

%%% integrate code below correctly!
if isa(S,'emptySet')
    % intersection with the empty set yields the empty set
    P_out = polytope.empty(dim(P));
    return
end

if isa(S,'fullspace')
    % R^n does not impose additional constraints
    P_out = polytope(P);
    return
end

% call function with lower precedence
if isa(S,'contSet') && S.precedence < P.precedence
    P_out = and_(S,P,type,varargin{:});
    return
end

% re-order such that first argument is polytope
[P,S] = findClassArg(P,S,'polytope');

% polytope
if isa(S,'polytope')
    P_out = aux_and_polytope(P,S);
    return
end

% zonotope bundle
if isa(S,'zonoBundle')
    if representsa_(P,'halfspace',1e-12)
        P_out = aux_and_zonoBundle_halfspace(P,S);
    elseif representsa_(P,'conHyperplane',1e-12)
        P_out = aux_and_zonoBundle_conHyperplane(P,S);
    else
        P_out = aux_and_polytope(P,polytope(S));
    end
    return
end

% constrained zonotope
if isa(S,'conZonotope')
    if representsa_(P,'conHyperplane',1e-12)
        P_out = aux_and_conZonotope_conHyperplane(P,S);
    else
        P_out = aux_and_conZonotope(P,S);
    end
    return
end

% intersection with numeric: if the point is contained, the intersection is
% that point, otherwise empty
if isnumeric(S) && size(S,2) == 1
    if contains_(P,S,'exact',eps)
        P_out = polytope(S);
    else
        P_out = polytope.empty(dim(P));
    end
    return
end

% all other cases
try
    % convert second object to polytope
    S = polytope(S);
catch ME
    % no conversion operation implemented
    throw(CORAerror('CORA:noops',P,S));
end
P_out = aux_and_polytope(P,S);

end


% Auxiliary functions -----------------------------------------------------

function P_out = aux_and_polytope(P,S)

% quick check: fully empty object -> fullspace
if representsa_(S,'fullspace',0)
    P_out = polytope(P);
    return
elseif representsa_(P,'fullspace',0)
    P_out = polytope(S);
    return
end

% ensure that both sets have the halfspace representation
constraints(P);
constraints(S);

% compute intersection
P_out = polytope([P.A_.val; S.A_.val], [P.b_.val; S.b_.val], ...
    [P.Ae_.val; S.Ae_.val], [P.be_.val; S.be_.val]);

% set properties
% intersection with empty set is empty
if (~isempty(P.emptySet.val) && P.emptySet.val) ...
        || (~isempty(S.emptySet.val) && S.emptySet.val)
    P_out.emptySet.val = true;
end

% intersection with bounded set yields a bounded set
if (~isempty(P.bounded.val) && P.bounded.val) ...
        || (~isempty(S.bounded.val) && S.bounded.val)
    P_out.bounded.val = true;
end

% intersection with a degenerate set yields a degenerate set
if (~isempty(P.fullDim.val) && ~P.fullDim.val) ...
        || (~isempty(S.fullDim.val) && ~S.fullDim.val)
    P_out.fullDim.val = false;
end

end

function zB_out = aux_and_zonoBundle_halfspace(P,zB)

% construct basis orthogonal to halfspace normal vector
B = gramSchmidt(P.A');

% compute enclosing interval in transformed space
Z_ = B' * zB.Z{1};
I_ = interval(Z_);

% consider upper bound applied by halfspace constraint c*x <= d
lb = infimum(I_);
ub = supremum(I_);

ub(1) = P.b/norm(P.A');

I_ = interval(lb,ub);

% backtransformation to orginal space
Z = B * zonotope(I_);

% intersection
zB_out = and_(zB,Z,'exact');

end

function zB_out = aux_and_zonoBundle_conHyperplane(P,zB)

% Part 1: intersection with the hyperplane

% construct basis orthogonal to halfspace normal vector
B = gramSchmidt(P.Ae(1,:)');

% compute enclosing interval in transformed space
Z_ = B' * zB.Z{1};
I_ = interval(Z_);

% consider upper bound applied by halfspace constraint c*x <= d
lb = infimum(I_);
ub = supremum(I_);

ub(1) = P.be(1)/norm(P.Ae(1,:)');
lb(1) = ub(1);

I_ = interval(lb,ub);

% backtransformation to orginal space
Z = B * zonotope(I_);

% intersection
zB_out = and_(zB,Z,'exact');


% Part 2: intersection with the constraints

% for each inequality constraint, check if the set is fully contained,
% and if not, intersect the set with that halfspace
C = P.A; d = P.b; nrCon = length(d);
for i = 1:nrCon
    P = polytope(C(i,:)',d(i));
    if ~contains_(P,zB_out,'exact',1e-8)
        zB_out = aux_and_zonoBundle_halfspace(P,zB_out);
    end
end

end

function cZ_out = aux_and_conZonotope_conHyperplane(P,cZ)

% calculate intersection between constrained zonotope and hyperplane
[P_A,P_b,P_Ae,P_be] = priv_normalizeConstraints(P.A_.val,P.b_.val,P.Ae_.val,P.be_.val,'A');
[P_A,P_b,P_Ae,P_be] = priv_compact_toEquality(P_A,P_b,P_Ae,P_be,1e-12);

c = cZ.c; G = cZ.G;
A = [cZ.A; P_Ae*G]; b = [cZ.b; P_be - P_Ae*c];
cZ_out = conZonotope([c,G],A,b); 

% for each inequality constraint, check if the set is fully contained,
% and if not, intersect the set with that halfspace
nrCon = numel(P_b);
for i = 1:nrCon
    hs = polytope(P_A(i,:),P_b(i));
    if ~contains_(hs,cZ_out,'exact',1e-6)
        cZ_out = aux_and_conZonotope(hs,cZ_out); 
    end
end

end
    
function cZ_out = aux_and_conZonotope(P,cZ)

G_ = cZ.G; c = cZ.c;
A_ = cZ.A; b_ = cZ.b;
n = dim(cZ);

% loop over inequality constraints
nrCon = length(P.b);
for i = 1:nrCon
    % current direction
    C = P.A(i,:); 

    % compute lower bound
    lb = supportFunc_(cZ,C','lower');
    % upper bound trivial
    ub = P.b(i);

    % extend matrices
    A_ = [A_, zeros(size(A_,1),1); C*G_, 0.5*(lb-ub)];
    b_ = [b_; 0.5*(ub+lb) - C*c];
    G_ = [G_, zeros(n,1)];
end

cZ_out = conZonotope(c,G_,A_,b_);

end

% ------------------------------ END OF CODE ------------------------------
