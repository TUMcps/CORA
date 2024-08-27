function res = isIntersecting_(P,S,type,varargin)
% isIntersecting_ - determines if a polytope intersects another set
%
% Syntax:
%    res = isIntersecting_(P,S,type)
%
% Inputs:
%    P - polytope object
%    S - contSet object, numerical vector, point cloud
%    type - type of check ('exact' or 'approx')
%    tol - tolerance
%
% Outputs:
%    res - true/false whether intersection occurs
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([-1 -1; -1 1; 1 1;0 -1],[0;2;2;0]) + [2;0];
%    P3 = P2 + [3;0];
%
%    isIntersecting(P1,P2)
%    isIntersecting(P1,P3)
%
%    figure; hold on
%    plot(P1,[1,2],'b');
%    plot(P2,[1,2],'g');
%
%    figure; hold on
%    plot(P1,[1,2],'b');
%    plot(P3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, interval/isIntersecting_

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       20-November-2019
% Last update:   14-December-2022 (MW, add call to MOSEK for linprogs)
%                27-April-2023 (MW, add support for points)
%                10-December-2023 (MW, fix for point clouds)
% Last revision: 10-July-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% set tolerance
tol = setDefaultValues({1e-12},varargin);

% get polytope object
[P,S] = findClassArg(P,S,'polytope');

% fullspace
if representsa_(P,'fullspace',0)
    % unless S is the empty set, there is an intersection
    res = ~representsa_(S,'emptySet',eps); return
elseif representsa_(S,'fullspace',0)
    res = ~representsa_(P,'emptySet',eps); return
end

% dimension
n = dim(P);

% 1D case: convert to intervals (fast) and check those
if n == 1
    try
        res = isIntersecting_(interval(P),interval(S),type,tol);
        return
    end
end


% nD cases
if isa(S,'polytope')
    % polytope and polytope
    res = aux_isIntersecting_poly_poly(P,S);
     
elseif strcmp(type,'exact')
    % exact check for intersection
        
    if isnumeric(S)
        % a point cloud intersects a polytope if any point is contained
        res = any(contains_(P,S,'exact',tol));
    elseif isa(S,'interval') || isa(S,'zonotope')
        res = aux_isIntersecting_P_cZ(P,conZonotope(S));
    elseif isa(S,'conZonotope')
        res = aux_isIntersecting_P_cZ(P,S);
    elseif isa(S,'zonoBundle')
        res = aux_isIntersecting_P_zB(P,S);
    elseif isa(S,'halfspace') || isa(S,'conHyperplane')
        res = isIntersecting_(S,P,'exact');
    else
        throw(CORAerror('CORAerror:noops',P,S));
    end
    
else
    % over-approximative check for intersection
    res = aux_isIntersecting_approx(P,S);
    
end

end
    

% Auxiliary functions -----------------------------------------------------

function res = aux_isIntersecting_poly_poly(P1,P2)

% go over possible cases in the following order:
% 1) H-polytope & H-polytope
% 2) H-polytope & V-polytope
% 3) V-polytope & V-polytope
if P1.isHRep.val
    if P2.isHRep.val
        % both in H
        res = aux_isIntersecting_Hpoly_HPoly(P1,P2);
    else
        % P1 in H, P2 in V
        res = aux_isIntersecting_Hpoly_Vpoly(P1,P2);
    end
elseif P2.isHRep.val
    % P1 in V, P2 in H (caution: input arguments switched!)
    res = aux_isIntersecting_Hpoly_Vpoly(P2,P1);
else
    % both in V
    res = aux_isIntersecting_Vpoly_Vpoly(P1,P2);
end


end

function res = aux_isIntersecting_Hpoly_HPoly(P1,P2)
% intersection check of two polytopes in halfspace representation

% construct the intersection
P1P2 = polytope([P1.A_.val; P2.A_.val], [P1.b_.val; P2.b_.val], ...
                [P1.Ae_.val; P2.Ae_.val], [P1.be_.val; P2.be_.val]);

% check if it is empty
res = ~representsa(P1P2,'emptySet');

end

function res = aux_isIntersecting_Hpoly_Vpoly(P1,P2)
% intersection check of a polytope in halfspace representation (P1) and
% another polytope in vertex representation (P2), see [1, (14)-(15)]

% min_{x,y} 1
% s.t.  A x <= b,
%       Ae x == be,
%       V beta = x,
%       sum_j beta_j = 1
%       beta_j >= 0

n = dim(P1);
nrIneq = length(P1.b_.val);
nrEq = length(P1.be_.val);
numVert = size(P2.V_.val,2);

problem.f = zeros(n+numVert,1);
problem.Aineq = [P1.A_.val, zeros(nrIneq,numVert);
                 zeros(numVert,n), -eye(numVert)];
problem.bineq = [P1.b_.val; zeros(numVert,1)];
problem.Aeq = [P1.Ae_.val, zeros(nrEq,numVert);
               zeros(1,n), ones(1,numVert)];
problem.beq = [P1.be_.val; 1];
problem.lb = [];
problem.ub = [];

% solve linear program
[~,val,exitflag] = CORAlinprog(problem);

% only relevant if it's infeasible or feasible
res = exitflag ~= -2;

end

function res = aux_isIntersecting_Vpoly_Vpoly(P1,P2)
% intersection check of two polytopes in vertex representation according to
% [1, (11)-(12)]

% min_{beta1, beta2} 1
% s.t.  V1 beta1 = V2 beta2
%       sum_j beta1_j = 1
%       beta1_j >= 0
%       sum_j beta2_j = 1
%       beta2_j >= 0

n = dim(P1);
numVert1 = size(P1.V,2);
numVert2 = size(P2.V,2);

problem.f = zeros(numVert1+numVert2,1);
problem.Aineq = [-eye(numVert1,numVert1), zeros(numVert1,numVert2);
                 zeros(numVert2,numVert1), -eye(numVert2)];
problem.bineq = [zeros(numVert1,1); zeros(numVert2,1)];
problem.Aeq = [P1.V, -P2.V;
               ones(1,numVert1), zeros(1,numVert2);
               zeros(1,numVert1), ones(1,numVert2)];
problem.beq = [zeros(n,1); 1; 1];
problem.lb = [];
problem.ub = [];

% solve linear program
[~,val,exitflag] = CORAlinprog(problem);

% only relevant if it's infeasible or feasible
res = exitflag ~= -2;

end

function res = aux_isIntersecting_P_cZ(P,cZ)
% check if a polytope {x | H*x <= d, He*x = de} and a constraint zonotope 
% {x = c + G*a | A*a = b, a \in [-1,1]} intersect by solving the following
% linear program:
%
% min_{t,x,a} t
%
% s.t. Hx - d <= t
%        He x = de
%      c + Ga = x
%          Aa = b
%           a \in [-1,1]

    % halfspaces and offsets of polytope
    H = P.A_.val;
    d = P.b_.val;
    He = P.Ae_.val;
    de = P.be_.val;
    
    % center and generator matrix of constrained zonotope
    c = cZ.c;
    G = cZ.G;
    
    % dimension
    n = length(c);
    % number of generators
    m = size(G,2);
    % number of polytope inequality constraints
    p = size(H,1);
    % number of polytope equality constraints
    pe = size(He,1);

    % optimization variable is [t;x;a] with length 1 + n + m

    % construct inequality constraints
    %   Hx - d <= t   =>  Hx <= t + d
    %   a \in [-1,1]  =>  -a <= 1, a <= 1
    
    % compute number of elements and minimum sparsity
    numElem = (p+2*m)*(1+n+m);
    minNumZeros = 2*m*(1+n) + p*m + 2*m*(m-1);
    sparsity = minNumZeros / numElem;
    % use sparse representation if beneficial
    if sparsity > 0.5
        Aineq = [-ones(p,1), H, sparse(p,m);
             sparse(m,1+n), -speye(m);
             sparse(m,1+n), speye(m)];
    else
        Aineq = [-ones(p,1), H, zeros(p,m);
             zeros(m,1+n), -eye(m);
             zeros(m,1+n), eye(m)];
    end
    bineq = [d;ones(2*m,1)];

    % construct equality constraints
    % He*x = de
    % x = c + G*a  <=>  x - G*a = c
    % A*a = b
    Aeq = [[zeros(pe,1),He,zeros(pe,m)]; [zeros(n,1),eye(n),-G]];
    beq = [de; c];
    
    if ~isempty(cZ.A)
        A = cZ.A; b = cZ.b;
        % A*a = b
        Aeq = [Aeq;
               zeros(size(A,1),1+n),A];
        beq = [beq;b];
    end
    
    % construct objective function
    % min t
    f = [1;zeros(n+m,1)];

    % init linprog struct
    problem.f = f;
    problem.Aineq = Aineq;
    problem.bineq = bineq;
    problem.Aeq = Aeq;
    problem.beq = beq;
    problem.lb = [];
    problem.ub = [];

    % solve linear program
    [~,val,exitflag] = CORAlinprog(problem);

    % multiple cases:
    % 1. constrained zonotope (and possibly polytope) is empty 
    %    => infeasible (dual unbounded)
    % 2. polytope empty => max(y)>0
    %     OR
    %    no intersection point between non-empty polytope and non-empty cZ
    % ELSE: error
    
    % cannot become unbounded above (left to the reader as an exercise);
    % also cannot become unbounded below: That implies there exists an
    % x in the intersection, and since cZ is bounded, so is H*x-d,
    % which is a lower bound on y

    if exitflag == -2
        % 1. no feasible point was found
        res = false;
    elseif exitflag == 1 % one could maybe include 3 and 0
        % 2. feasible point found, check if max(y) < 0 == val
        % tol? => use constraint tolerance (default: 1e-6)
        tol = 1e-6;
        res = val < tol;
    else
        throw(CORAerror('CORA:solverIssue','linprog'));
    end

end
    
function res = aux_isIntersecting_P_zB(P,zB)
% check if a polytope {x | H*x <= d} and a zonotope bundle 
% {x = c1 + G1*a|a \in [-1,1]} \cup ... \cup {x = cq + Gq*a | a \in [-1,1]}
% intersect by solving the following linear program:
%
% min sum(y)
%
% s.t. Hx - y < d
%           y > 0
%  c1 + G1 a1 = x
%             .
%             .
%             .
%  cq + Gq a2 = x
%   a1,...,aq \in [-1,1]

    % get object properties
    H = P.A_.val;
    d = P.b_.val;
    
    [p,n] = size(H);

    % construct inequality constraints
    A = [H,-eye(p);zeros(p,n),-eye(p)];
    b = [d;zeros(p,1)];
    
    % loop over all parallel zonotopes in the bundle
    Aeq = [];
    beq = [];
    
    for i = 1:zB.parallelSets
       
        Z = zB.Z{i};
        c = center(Z);
        G = generators(Z);
        m = size(G,2);
        
        % construct equality constraints 
        Aeq = blkdiag(Aeq,-G);
        beq = [beq;c];
        
        % construct inequality constraints
        A = blkdiag(A,[eye(m);-eye(m)]);
        b = [b;ones(2*m,1)];
    end
    
    temp = repmat(eye(n),[zB.parallelSets,1]);
    Aeq = [temp,zeros(size(temp,1),p),Aeq];
    
    % construct objective function
    f = zeros(size(Aeq,2),1);
    f(n+1:n+p) = ones(p,1);

    % init linprog struct
    problem.f = f;
    problem.Aineq = Aineq;
    problem.bineq = bineq;
    problem.Aeq = Aeq;
    problem.beq = beq;
    problem.lb = [];
    problem.ub = [];

    % solve linear program
    [~,val,exitflag] = CORAlinprog(problem);

    % check if intersection between the two sets is empty
    if exitflag < 0 || (val > 0 && ~withinTol(val,0))
        res = false; 
    else
        res = true;
    end
    
end

function res = aux_isIntersecting_approx(P,S)
% approx check, i.e., conservative intersection check (possibility of false
% positives)
        
res = true;
A = P.A_.val;
b = P.b_.val;

% special 'approx' algorithm for zonotope bundles
if isa(S,'zonoBundle')

    % loop over all parallel zonotopes
    for j = 1:S.parallelSets
        
        % read out j-th zonotope
        Z = S.Z{j};
        
        % loop over all halfspaces
        for i = 1:size(A,1)
            if b(i) < supportFunc_(Z,A(i,:)','lower')
                res = false; return;
            end
        end

    end
    
else

    otherOptions = {};
    if isa(S,'conPolyZono') || isa(S,'polyZonotope')
        otherOptions = {'interval',8,1e-3};
    end
    
    % loop over all halfspaces
    for i = 1:size(A,1)
        if b(i) < supportFunc_(S,A(i,:)','lower',otherOptions{:})
            res = false; return;
        end
    end

end
end

% ------------------------------ END OF CODE ------------------------------
