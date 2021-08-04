function res = isIntersecting(obj1,obj2,varargin)
% isIntersecting - determines if constrained hyperplane obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%    res = isIntersecting(obj1,obj2,type)
%
% Inputs:
%    obj1 - conHyperplane object
%    obj2 - contSet object
%    type - type of check ('exact' or 'approx')
%
% Example: 
%    hp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    zono = zonotope([0 1 1 0; 0 1 0 1]);
%    poly = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]) + [2;2];
%
%    isIntersecting(hp,zono)
%    isIntersecting(hp,poly)
% 
%    figure
%    hold on
%    xlim([-6,6]);
%    ylim([-6,6]);
%    plot(hp,[1,2],'b');
%    plot(zono,[1,2],'g');
% 
%    figure
%    hold on
%    xlim([-6,6]);
%    ylim([-6,6]);
%    plot(hp,[1,2],'b');
%    plot(poly,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/isIntersecing

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  14-Sep-2019
%               20-Nov-2019
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && ~isempty(varargin{1}) 
        type = varargin{1};
    end
    
    % get conHyperplane object
    if ~isa(obj1,'conHyperplane')
       temp = obj1;
       obj1 = obj2;
       obj2 = temp;
    end
    
    % exact or apprxomiative algorithm
    if strcmp(type,'exact')
        
        if isempty(obj1.C) && ~isa(obj2,'zonoBundle')
            res = isIntersecting(obj1,obj2,'approx');
        else
            if isa(obj2,'mptPolytope')
                res = intersectPolyPoly(obj1,obj2);
            elseif isa(obj2,'interval') || isa(obj2,'zonotope')
                res = intersectPolyConZono(obj1,conZonotope(obj2));
            elseif isa(obj2,'conZonotope')
                res = intersectPolyConZono(obj1,obj2);
            elseif isa(obj2,'zonoBundle')
                res = intersectPolyZonoBundle(obj1,obj2);  
            else
                error('No exact algorithm implemented for this set representation!');
            end
        end
        
    else
       
        res = 1;
        
        % special 'approx' algorithm for zonotope bundles
        if isa(obj2,'zonoBundle')
        
            % loop over all parallel zonotopes
            for j = 1:length(obj2.Z)

                zono = zonotope(obj2.Z{j});
            
                % check instesection with hyperplane
                infi = supportFunc(zono,obj1.h.c,'lower');
                sup = supportFunc(zono,obj1.h.c,'upper');

                if ~in(interval(infi,sup),obj1.h.d)
                   res = 0;
                   return;
                end

                % check intersection with inequality constraints
                C = obj1.C;
                d = obj1.d;

                for i = 1:size(C,1)
                   bound = supportFunc(zono,C(i,:)','lower');
                   if bound > d(i)
                      res = 0;
                      return;
                   end
                end
            end 
            
        else
            
            % check instesection with hyperplane
            infi = supportFunc(obj2,obj1.h.c,'lower');
            sup = supportFunc(obj2,obj1.h.c,'upper');

            if ~in(interval(infi,sup),obj1.h.d)
               res = 0;
               return;
            end

            % check intersection with inequality constraints
            C = obj1.C;
            d = obj1.d;

            for i = 1:size(C,1)
               bound = supportFunc(obj2,C(i,:)','lower');
               if bound > d(i)
                  res = 0;
                  return;
               end
            end
        end
    end
end


% Auxiliary Functions -----------------------------------------------------    

function res = intersectPolyPoly(obj1,obj2)
% check if a constrained hyperplane {x | c x = d, A x < b} and a polytope
% {x | H x <= k} intersect by solving the following linear program
%
% min sum(y)
%
% s.t. H x - y < k
%      A x - y < b
%            y > 0
%          c x = d 
    
    % construct matrices for inequality constraints
    A = [obj2.P.A;obj1.C];
    b = [obj2.P.b;obj1.d];
    
    % construct matrices for equality constraint
    Aeq = obj1.h.c';
    beq = obj1.h.d;
    
    % introduce slack variables y
    m = length(b);
    n = size(A,2);
    
    A = [A,-eye(m);zeros(m,n),-eye(m)];
    b = [b;zeros(m,1)];

    Aeq = [Aeq,zeros(1,m)];
    
    f = [zeros(n,1);ones(m,1)];
    
    % solve the dual problem using linear programming
    options = optimoptions('linprog','display','off');
    
    [~,val,exitflag] = linprog(f,A,b,Aeq,beq,[],[],options);

    % check if intersection between the two polytopes is empty
    res = 1;
    
    if exitflag < 0 || val > eps
       res = 0; 
    end
end

function res = intersectPolyConZono(obj1,obj2)
% check if a constrained hyperplane {x | h x = d, H x < k}   
% and a constraint zonotope {x = c + G*a | A*a = b, a \in [-1,1]} intersect 
% by solving the following linear program:
%
% min sum(y)
%
% s.t.     cx = d
%      Hx - y < k
%           y > 0
%      c + Ga = x
%          Aa = b
%           a \in [-1,1]

    % get object properties
    h = obj1.h.c;
    d = obj1.h.d;
    
    H = obj1.C;
    k = obj1.d;
    
    c = obj2.Z(:,1);
    G = obj2.Z(:,2:end);
    
    n = length(c);
    m = size(G,2);

    % construct inequality constraints
    if isempty(H)
       H = [h';-h'];
       k = [d;-d];
    else
       H = [H;h';-h'];
       k = [k;d;-d]; 
    end 
    
    p = size(H,1);
    
    A = [H,-eye(p);zeros(p,n),-eye(p)];
    b = [k;zeros(p,1)];
    
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
% check if a constrained hyperplane {x | h x = d, H x < k} and a bundle
% {x = c1 + G1*a|a \in [-1,1]} \cup ... \cup {x = cq + Gq*a | a \in [-1,1]}
% intersect by solving the following linear program:
%
% min sum(y)
%
% s.t.     hx = d
%      Hx - y < k
%           y > 0
%  c1 + G1 a1 = x
%             .
%             .
%             .
%  cq + Gq a2 = x
%   a1,...,aq \in [-1,1]

    % get object properties
    h = obj1.h.c;
    d = obj1.h.d;
    
    H = obj1.C;
    k = obj1.d;

    % construct inequality constraints
    if isempty(H)
       H = [h';-h'];
       k = [d;-d];
    else
       H = [H;h';-h'];
       k = [k;d;-d]; 
    end 
    
    [p,n] = size(H);
    
    A = [H,-eye(p);zeros(p,n),-eye(p)];
    b = [k;zeros(p,1)];
    
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