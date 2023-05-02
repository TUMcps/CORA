function res = isIntersecting_(P,S,type,varargin)
% isIntersecting_ - determines if a polytope intersects a set
%
% Syntax:  
%    res = isIntersecting_(P,S)
%    res = isIntersecting_(P,S,type)
%
% Inputs:
%    P - mptPolytope object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    P1 = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = mptPolytope([-1 -1; -1 1; 1 1;0 -1],[0;2;2;0]) + [2;0];
%    P3 = P2 + [3;0];
%
%    isIntersecting(P1,P2)
%    isIntersecting(P1,P3)
%
%    figure; hold on;
%    plot(P1,[1,2],'b');
%    plot(P2,[1,2],'g');
%
%    figure; hold on;
%    plot(P1,[1,2],'b');
%    plot(P3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/isIntersecting_

% Author:       Niklas Kochdumper, Victor Gassmann
% Written:      20-November-2019
% Last update:  27-April-2023 (VG: changed LPs to more robust version)
% Last revision:27-March-2023 (MW, rename isIntersecting_)

%------------- BEGIN CODE --------------
   
    % init options for linear programs
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end

    % polytope and polytope
    if isa(S,'mptPolytope')
        
        res = aux_intersectPolyPoly(P,S,options);
         
    else
        
        % exact check for intersection
        if strcmp(type,'exact')
            
            if isa(S,'interval') || isa(S,'zonotope')
                res = aux_intersectPolyConZono(P,conZonotope(S),options);
            elseif isa(S,'conZonotope')
                res = aux_intersectPolyConZono(P,S,options);
            elseif isa(S,'zonoBundle')
                res = aux_intersectPolyZonoBundle(P,S,options);
            elseif isa(S,'halfspace') || isa(S,'conHyperplane')
                res = isIntersecting_(S,P,type);
            else
                throw(CORAerror('CORA:noops',P,S));
            end
            
        % over-approximative check for intersection    
        else
            
            res = true;
            A = P.P.A;
            b = P.P.b;
            
            % special 'approx' algorithm for zonotope bundles
            if isa(S,'zonoBundle')

                % loop over all parallel zonotopes
                for j = 1:length(S.Z)
                    
                    Z = S.Z{j};
                    
                    % loop over all halfspaces
                    for i = 1:size(P.P.A,1)

                        val = supportFunc_(Z,A(i,:)','lower');

                        if val > b(i)
                           res = false;
                           return;
                        end
                    end
                end
                
            else

                otherOptions = {};
                if isa(S,'polyZonotope') || isa(S,'conPolyZono')
                    otherOptions = {'interval',8,1e-3};
                end
                
                % loop over all halfspaces
                for i = 1:size(P.P.A,1)

                    val = supportFunc_(S,A(i,:)','lower',otherOptions{:});

                    if val > b(i)
                       res = false;
                       return;
                    end
                end
            end
        end
    end
end
    

% Auxiliary Functions -----------------------------------------------------    

function res = aux_intersectPolyPoly(P1,P2,options)
% check if two polytopes intersect
    
    % construct matrices for inequality constraints A * x <= b
    A = [P1.P.A;P2.P.A];
    b = [P1.P.b;P2.P.b];
    
    % introduce slack variables y to solve the dual problem
    % min t
    % s.t. Ax - b <= y,
    %         y_i <= t
    
    % opt var x = [t;y;x]

    p = length(b);
    n = size(A,2);

    % ineqs
    % Ax -b <= y
    % y_i <= t
    Aineq = [zeros(p,1),-eye(p),A;
             -ones(p,1),eye(p),zeros(p,n)];
    bineq = [b;zeros(p,1)];
    
    f = [1;zeros(p+n,1)];
    
    % solve the dual problem using linear programming
    [~,val,exitflag] = linprog(f,Aineq,bineq,[],[],[],[],options);

    % polytope(s) empty => max(y)>0 or there is the possibilty of
    % intersection
    % ELSE: error
    
    % cannot become unbounded above (obviously)

    if exitflag==1 % one could maybe include 3 and 0
        % 2.
        % check if max(y) < 0 == val
        % tol? => use constraint tolerance
        tol = options.ConstraintTolerance;
        res = val < tol;
    elseif exitflag==-3
        res = true;
    else
        throw(CORAerror('CORA:solverIssue','linprog'));
    end
end

function res = aux_intersectPolyConZono(P,cZ,options)
% check if a polytope {x | H*x <= d} and a constraint zonotope 
% {x = c + G*a | A*a = b, a \in [-1,1]} intersect by solving the following
% linear program:
%
% min t
%
% s.t. Hx -d <= y
%        y_i <= t, 1\leq i\leq p (number of rows in H),
%      c + Ga = x
%          Aa = b
%           a \in [-1,1]

    % get object properties
    H = P.P.A;
    d = P.P.b;
    
    c = cZ.Z(:,1);
    G = cZ.Z(:,2:end);
    
    n = length(c);
    m = size(G,2);
    p = size(H,1);

    % opt var x = [t;y;x;a] 

    % construct inequality constraints
    % Hx -d <=y
    % y_i <= t
    % -a <= 1
    %  a <= 1
    Aineq = [zeros(p,1),-eye(p),H,zeros(p,m);
             -ones(p,1),eye(p),zeros(p,n+m);
             zeros(m,1+p+n),-eye(m);
             zeros(m,1+p+n),eye(m)];
    bineq = [d;zeros(p,1);ones(2*m,1)];

    % construct equality constraints
    % x = c+ G*a
    % A*a = b
    Aeq = [zeros(n,1+p),eye(n),-G];
    beq = c;
    
    if ~isempty(cZ.A)
        A = cZ.A;
        b = cZ.b;
        % A*a = b
        Aeq = [Aeq;
               zeros(size(A,1),1+p+n),A];
        beq = [beq;b];
    end
    
    % construct objective function
    % min t
    f = [1;zeros(p+n+m,1)];
    
    % solve linear program
    [~,val,exitflag] = linprog(f,Aineq,bineq,Aeq,beq,[],[],options);

    % multiple cases:
    % 1. constrained zonotope (and possibly polytope) is empty 
    %    => infeasible (dual unbounded)
    % 2. polytope empty => max(y)>0
    %     OR
    %    no intersection point between non-empty polytope and non-empty cZ
    % ELSE: error
    
    % cannot become unbounded above (obviously); also cannot become
    % unbounded below: That implies there exists an x in the intersection,
    % and since cZ is bounded, so is H*x-d, which is a lower bound on y

    if exitflag==-2
        % 1.
        res = false;
    elseif exitflag==1 % one could maybe include 3 and 0
        % 2.
        % check if max(y) < 0 == val
        % tol? => use constraint tolerance
        tol = options.ConstraintTolerance;
        res = val < tol;
    else
        throw(CORAerror('CORA:solverIssue','linprog'));
    end

end
    
function res = aux_intersectPolyZonoBundle(P,zB,options)
% check if a polytope {x | H*x <= d} and a zonotope bundle 
% {x = c1 + G1*a|a \in [-1,1]} \cap ... \cap {x = cq + Gq*a | a \in [-1,1]}
% intersect by solving the following linear program:
%
% min t
%
% s.t. Hx - d <= y,
%         y_i <= t,
%  c1 + G1 a1 = x,
%             .
%             .
%             .
%  cq + Gq a2 = x,
%   a1,...,aq \in [-1,1]

    % get object properties
    H = P.P.A;
    d = P.P.b;
    
    n = dim(P);
    p = size(H,1);

    % opt var x = [t;y;x;aa]
    

    % construct inequality constraints
    % Hx -d <=y
    % y_i <= t
    % do not include aa yet (don't know how large aa is yet)
    Aineq_ = [zeros(p,1),-eye(p),H;
             -ones(p,1),eye(p),zeros(p,n)];
    bineq = [d;zeros(p,1)];
    
    mm = 0;
    N = zB.parallelSets;
    beq = zeros(N*n,1);
    
    GG_c = cell(N,1);
    for i = 1:N
       
        Z = zB.Z{i};
        c = center(Z);
        G = generators(Z);
        m = size(G,2);
        mm = mm +m;
        GG_c{i} = G;
        beq((i-1)*n+(1:n)) = c;
    end
    
    % extend ineq constraints by length of aa
    Aineq = [Aineq_,zeros(size(Aineq_,1),mm)];
   
    % construct equality constraints
    % x = c_k +G_k*a_k
    Aeq = [zeros(N*n,1+p),repmat(eye(n),N,1),-blkdiag(GG_c{:})];

    % lower and upper bounds
    lb = [-inf(1+p+n,1);-ones(mm,1)];
    ub = [inf(1+p+n,1);ones(mm,1)];
    
    % construct objective function
    % min t
    f = [1;zeros(p+n+mm,1)];
    
    % solve linear program
    [~,val,exitflag] = linprog(f,Aineq,bineq,Aeq,beq,lb,ub,options);
    
    % multiple cases:
    % 1. zonoBundle (and possibly polytope) is empty 
    %    => infeasible (dual unbounded)
    % 2. polytope empty => max(y)>0
    %     OR
    %    no intersection point between non-empty polytope and non-empty cZ
    % ELSE: error
    
    % cannot become unbounded above (obviously); also cannot become
    % unbounded below: That implies there exists an x in the intersection,
    % and since cZ is bounded, so is H*x-d, which is a lower bound on y

    if exitflag==-2
        % 1.
        res = false;
    elseif exitflag==1 % one could maybe include 3 and 0
        % 2.
        % check if max(y) < 0 == val
        % tol? => use constraint tolerance
        tol = options.ConstraintTolerance;
        res = val < tol;
    else
        throw(CORAerror('CORA:solverIssue','linprog'));
    end
end
    
%------------- END OF CODE --------------