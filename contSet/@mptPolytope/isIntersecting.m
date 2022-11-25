function res = isIntersecting(obj1,obj2,varargin)
% isIntersecting - determines if a the sets obj1 and obj2 intersect
%
% Syntax:  
%    res = isIntersecting(obj1,obj2,type)
%    res = isIntersecting(obj1,obj2)
%
% Inputs:
%    obj1 - contSet object
%    obj2 - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    poly1 = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    poly2 = mptPolytope([-1 -1; -1 1; 1 1;0 -1],[0;2;2;0]) + [2;0];
%    poly3 = poly2 + [3;0];
%
%    isIntersecting(poly1,poly2)
%    isIntersecting(poly1,poly3)
%
%    figure
%    hold on
%    plot(poly1,[1,2],'b');
%    plot(poly2,[1,2],'g');
%
%    figure
%    hold on
%    plot(poly1,[1,2],'b');
%    plot(poly3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/isIntersecting

% Author:       Niklas Kochdumper
% Written:      20-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && ~isempty(varargin{1}) 
        type = varargin{1};
    end
    
    % get mptPolytope object
    if ~isa(obj1,'mptPolytope')
       temp = obj1;
       obj1 = obj2;
       obj2 = temp;
    end
    
    % polytope and polytope
    if isa(obj2,'mptPolytope')
        
        res = intersectPolyPoly(obj1,obj2);
         
    else
        
        % exact check for intersection
        if strcmp(type,'exact')
            
            if isa(obj2,'interval') || isa(obj2,'zonotope')
                res = intersectPolyConZono(obj1,conZonotope(obj2));
            elseif isa(obj2,'conZonotope')
                res = intersectPolyConZono(obj1,obj2);
            elseif isa(obj2,'zonoBundle')
                res = intersectPolyZonoBundle(obj1,obj2);  
            elseif isa(obj2,'halfspace') || isa(obj2,'conHyperplane')
                res = isIntersecting(obj2,obj1);
            else
                error('No exact algorithm implemented for this set representation!');
            end
            
        % over-approximative check for intersection    
        else
            
            res = 1;
            A = obj1.P.A;
            b = obj1.P.b;
            
            % special 'approx' algorithm for zonotope bundles
            if isa(obj2,'zonoBundle')

                % loop over all parallel zonotopes
                for j = 1:length(obj2.Z)
                    
                    zono = zonotope(obj2.Z{j});
                    
                    % loop over all halfspaces
                    for i = 1:size(obj1.P.A,1)

                        val = supportFunc(zono,A(i,:)','lower');

                        if val > b(i)
                           res = 0;
                           return;
                        end
                    end
                end
                
            else
                
                % loop over all halfspaces
                for i = 1:size(obj1.P.A,1)

                    val = supportFunc(obj2,A(i,:)','lower');

                    if val > b(i)
                       res = 0;
                       return;
                    end
                end
            end
        end
    end
end
    

% Auxiliary Functions -----------------------------------------------------    

function res = intersectPolyPoly(obj1,obj2)
% check if two polytopes intersect
    
    % construct matrices for inequality constraints A * x <= b
    A = [obj1.P.A;obj2.P.A];
    b = [obj1.P.b;obj2.P.b];
    
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
    options = optimoptions('linprog','display','off');
    
    [~,val,exitflag] = linprog(f,A,b,[],[],[],[],options);

    % check if intersection between the two polytopes is empty
    res = 1;
    
    if exitflag < 0 || val > eps
       res = 0; 
    end
end

function res = intersectPolyConZono(obj1,obj2)
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
    H = obj1.P.A;
    d = obj1.P.b;
    
    c = obj2.Z(:,1);
    G = obj2.Z(:,2:end);
    
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
    
    if ~isempty(obj2.A)
        Aeq = [Aeq;[zeros(size(obj2.A,1),n+p),obj2.A]];
        beq = [beq;obj2.b];
    end
    
    % construct objective function
    f = [zeros(n,1);ones(p,1);zeros(m,1)];
    
    % solve linear program
    options = optimoptions('linprog','display','off');
    
    [~,val,exitflag] = linprog(f,A,b,Aeq,beq,[],[],options);

    % check if intersection between the two polytopes is empty
    res = 1;
    
    if exitflag < 0 || val > eps
       res = 0; 
    end
end
    
function res = intersectPolyZonoBundle(obj1,obj2)
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
    H = obj1.P.A;
    d = obj1.P.b;
    
    [p,n] = size(H);

    % construct inequality constraints
    A = [H,-eye(p);zeros(p,n),-eye(p)];
    b = [d;zeros(p,1)];
    
    % loop over all parallel zonotopes in the bundle
    Aeq = [];
    beq = [];
    
    for i = 1:obj2.parallelSets
       
        zono = obj2.Z{i};
        c = center(zono);
        G = generators(zono);
        m = size(G,2);
        
        % construct equality constraints 
        Aeq = blkdiag(Aeq,-G);
        beq = [beq;c];
        
        % construct inequality constraints
        A = blkdiag(A,[eye(m);-eye(m)]);
        b = [b;ones(2*m,1)];
    end
    
    temp = repmat(eye(n),[obj2.parallelSets,1]);
    Aeq = [temp,zeros(size(temp,1),p),Aeq];
    
    % construct objective function
    f = zeros(size(Aeq,2),1);
    f(n+1:n+p) = ones(p,1);
    
    % solve linear program
    options = optimoptions('linprog','display','off');
    
    [~,val,exitflag] = linprog(f,A,b,Aeq,beq,[],[],options);

    % check if intersection between the two polytopes is empty
    res = 1;
    
    if exitflag < 0 || val > eps
       res = 0; 
    end
end


    
%------------- END OF CODE --------------