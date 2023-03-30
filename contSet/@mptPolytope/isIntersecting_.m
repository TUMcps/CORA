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

% Author:       Niklas Kochdumper
% Written:      20-November-2019
% Last update:  ---
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
    % min sum(y)
    % s.t. Ax - y <= b
    %           y >= 0
    
    m = length(b);
    n = size(A,2);
    A = [A,-eye(m);zeros(m,n),-eye(m)];
    b = [b;zeros(m,1)];
    
    f = [zeros(n,1);ones(m,1)];
    
    % solve the dual problem using linear programming
    [~,val,exitflag] = linprog(f,A,b,[],[],[],[],options);

    % check if intersection between the two polytopes is empty
    res = true;
    
    if exitflag < 0 || val > eps
       res = false; 
    end
end

function res = aux_intersectPolyConZono(P,cZ,options)
% check if a polytope {x | H*x <= d} and a constraint zonotope 
% {x = c + G*a | A*a = b, a \in [-1,1]} intersect by solving the following
% linear program:
%
% min sum(y)
%
% s.t. Hx - y < d
%           y > 0
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

    % construct inequality constraints
    A = [H,-eye(p);zeros(p,n),-eye(p)];
    b = [d;zeros(p,1)];
    
    A = blkdiag(A,[eye(m);-eye(m)]);
    b = [b;ones(m,1);ones(m,1)];

    % construct equality constraints
    Aeq = [eye(n),zeros(n,p),-G];
    beq = c;
    
    if ~isempty(cZ.A)
        Aeq = [Aeq;[zeros(size(cZ.A,1),n+p),cZ.A]];
        beq = [beq;cZ.b];
    end
    
    % construct objective function
    f = [zeros(n,1);ones(p,1);zeros(m,1)];
    
    % solve linear program
    [~,val,exitflag] = linprog(f,A,b,Aeq,beq,[],[],options);

    % check if intersection between the two polytopes is empty
    res = true;
    
    if exitflag < 0 || val > eps
        res = false; 
    end
end
    
function res = aux_intersectPolyZonoBundle(P,zB,options)
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
    H = P.P.A;
    d = P.P.b;
    
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
    
    % solve linear program
    [~,val,exitflag] = linprog(f,A,b,Aeq,beq,[],[],options);

    % check if intersection between the two polytopes is empty
    res = true;
    
    if exitflag < 0 || val > eps
        res = false; 
    end
end
    
%------------- END OF CODE --------------