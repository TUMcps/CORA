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
% See also: contSet/isIntersecting, interval/isIntersecting

% Authors:       Niklas Kochdumper
% Written:       20-November-2019
% Last update:   14-December-2022 (MW, add call to MOSEK for linprogs)
%                27-April-2023 (MW, add support for points)
%                10-December-2023 (MW, fix for point clouds)
% Last revision: ---

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

persistent isMosek
if isempty(isMosek)
    isMosek = isSolverInstalled('mosek');
end

persistent options
if isempty(options)
    options = optimoptions('linprog','display','off');
end

% dimension
n = dim(P);

% 1D case: convert to interval and check there
if n == 1
    res = isIntersecting_(interval(P),interval(S),type,tol);
    return
end


% polytope and polytope
if isa(S,'polytope')
    
    res = aux_intersectPolyPoly(P,S);
    %res = aux_intersectPolyPoly(P,S,isMosek,options);
     
else
    
    % exact check for intersection
    if strcmp(type,'exact')
        
        if isnumeric(S)
            % a point cloud intersects a polytope if any point is contained
            res = any(contains_(P,S,'exact',tol));
        elseif isa(S,'interval') || isa(S,'zonotope')
            res = aux_intersectPolyConZono(P,conZonotope(S),isMosek,options);
        elseif isa(S,'conZonotope')
            res = aux_intersectPolyConZono(P,S,isMosek,options);
        elseif isa(S,'zonoBundle')
            res = aux_intersectPolyZonoBundle(P,S,isMosek,options);
        elseif isa(S,'halfspace') || isa(S,'conHyperplane')
            res = isIntersecting_(S,P,'exact');
        else
            throw(CORAerror('CORAerror:noops',P,S));
        end
        
    % over-approximative check for intersection    
    else
        
        res = true;
        A = P.A;
        b = P.b;
        
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
end

end
    

% Auxiliary functions -----------------------------------------------------

function res = aux_intersectPolyPoly(P1,P2)

    % check dimensions
    if dim(P1) ~= dim(P2)
        error('The two polytopes need to be of the same dimensionality');
    end

    % construct the intersection
    res = true;
    P1P2 = polytope([P1.A; P2.A], [P1.b; P2.b], [P1.Ae; P2.Ae], [P1.be; P2.be]);

    % check if it is empty
    if representsa(P1P2, 'emptySet')
        res = false;
    end

end

function res = aux_intersectPolyConZono(P,cZ,isMosek,options)
% check if a polytope {x | H*x <= d} and a constraint zonotope 
% {x = c + G*a | A*a = b, a \in [-1,1]} intersect by solving the following
% linear program:
%
% min t
%
% s.t. Hx - d <= t
%      c + Ga = x
%          Aa = b
%           a \in [-1,1]

    % halfspaces and offset of polytope
    H = P.A;
    d = P.b;
    
    % center and generator matrix of constrained zonotope
    c = cZ.c;
    G = cZ.G;
    
    % dimension
    n = length(c);
    % number of generators
    m = size(G,2);
    % number of polytope halfspaces
    p = size(H,1);

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
    % x = c + G*a  =>  x - G*a = c
    % A*a = b
    Aeq = [zeros(n,1),eye(n),-G];
    beq = c;
    
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

    if isMosek

        % rewrite to MOSEK syntax
        a = [Aineq; Aeq];
        blc = [-Inf(p+2*m,1);beq];
        buc = [bineq;beq];

        resMOSEK = msklpopt(f,a,blc,buc,[],[],[],'minimize echo(0)');

        % check if intersection between the two sets is empty
        if strcmp(resMOSEK.sol.itr.prosta,'PRIMAL_INFEASIBLE')
            res = false;
        elseif strcmp(resMOSEK.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
            if resMOSEK.sol.itr.dobjval > 0 && ~withinTol(resMOSEK.sol.itr.dobjval,0)
                res = false;
            else
                res = true;
            end
        elseif strcmp(resMOSEK.sol.itr.prosta,'DUAL_INFEASIBLE')
            % dual cannot be infeasible, since problem cannot be unbounded
            throw(CORAerror('CORA:solverIssue',...
                'Problem is unbounded, but cannot be unbounded.'));
        end

    else
    
        % solve linear program
        [~,val,exitflag] = linprog(f,Aineq,bineq,Aeq,beq,[],[],options);

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
            % tol? => use constraint tolerance
            tol = options.ConstraintTolerance;
            res = val < tol;
        else
            throw(CORAerror('CORA:solverIssue','linprog'));
        end

    end

end
    
function res = aux_intersectPolyZonoBundle(P,zB,isMosek,options)
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
    H = P.A;
    d = P.b;
    
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
    
    if isMosek

        % rewrite to MOSEK syntax
        a = [A; Aeq];
        blc = [-Inf(size(A,1),1);beq];
        buc = [b;beq];
        
        % call MOSEK
        resMOSEK = msklpopt(f,a,blc,buc,[],[],[],'minimize echo(0)');

        % check if intersection is empty
        res = aux_checkIntersectionEmptyMOSEK(resMOSEK);

    else
        % MATLAB linprog

        % solve linear program
        [~,val,exitflag] = linprog(f,A,b,Aeq,beq,[],[],options);
        res = aux_checkIntersectionEmptyMATLAB(val,exitflag);
    end
    
end

function res = aux_checkIntersectionEmptyMOSEK(resMOSEK)
% same check for all helper functions

    % check if intersection between the two sets is empty
    if strcmp(resMOSEK.sol.itr.prosta,'PRIMAL_INFEASIBLE') ...
        || strcmp(resMOSEK.sol.itr.prosta,'DUAL_INFEASIBLE')
        res = false;
    elseif strcmp(resMOSEK.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
        if resMOSEK.sol.itr.dobjval > 0 || withinTol(resMOSEK.sol.itr.dobjval,0)
            res = false;
        else
            res = true;
        end
    end

end

function res = aux_checkIntersectionEmptyMATLAB(val,exitflag)
% same check for all helper functions

    % check if intersection between the two sets is empty
    if exitflag < 0 || (val > 0 && ~withinTol(val,0))
        res = false; 
    else
        res = true;
    end

end
    
% ------------------------------ END OF CODE ------------------------------
