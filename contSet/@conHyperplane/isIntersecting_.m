function res = isIntersecting_(hyp,S,type,varargin)
% isIntersecting_ - determines if a constrained hyperplane intersects a set
%
% Syntax:
%    res = isIntersecting_(hyp,S)
%    res = isIntersecting_(hyp,S,type)
%
% Inputs:
%    hyp - conHyperplane object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    hyp = conHyperplane([1 1],0,[1 0;-1 0],[2;2]);
%    Z = zonotope([0 1 1 0; 0 1 0 1]);
%    P = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]) + [2;2];
%
%    isIntersecting(hyp,Z)
%    isIntersecting(hyp,P)
% 
%    figure; hold on; xlim([-6,6]); ylim([-6,6]);
%    plot(hyp,[1,2],'b');
%    plot(Z,[1,2],'g');
% 
%    figure; hold on; xlim([-6,6]); ylim([-6,6]);
%    plot(hyp,[1,2],'b');
%    plot(P,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, polytope/isIntersecting_, conHyperplane/contains_

% Authors:       Niklas Kochdumper
% Written:       16-May-2018
% Last update:   14-September-2019
%                20-November-2019
%                18-July-2024 (MW, call polytope functions)
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------
    
    % exact or apprxomiative algorithm
    if strcmp(type,'exact')
        
        if isempty(hyp.C) && ~isa(S,'zonoBundle')
            res = isIntersecting(hyp,S,'approx');
        else

            if isa(S,'polytope')
                res = aux_intersectPolyPoly(hyp,S);
            elseif isa(S,'interval') || isa(S,'zonotope')
                res = isIntersecting_(polytope(hyp),conZonotope(S),type);
            elseif isa(S,'conZonotope')
                res = isIntersecting_(polytope(hyp),S,type);
            elseif isa(S,'zonoBundle')
                res = aux_intersectPolyZonoBundle(hyp,S);
            else
                throw(CORAerror('CORA:noExactAlg',S,type));
            end
        end
        
    else
       
        res = true;
        
        % special 'approx' algorithm for zonotope bundles
        if isa(S,'zonoBundle')
        
            % loop over all parallel zonotopes
            for j = 1:length(S.Z)

                Z = S.Z{j};
            
                % check instesection with hyperplane
                lb = supportFunc_(Z,hyp.a','lower');
                ub = supportFunc_(Z,hyp.a','upper');

                if ~contains_(interval(lb,ub),hyp.b,'exact',0)
                   res = false;
                   return;
                end

                % check intersection with inequality constraints
                C = hyp.C; d = hyp.d;

                for i = 1:size(C,1)
                   bound = supportFunc_(Z,C(i,:)','lower');
                   if bound > d(i)
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
            
            % check instesection with hyperplane
            lb = supportFunc_(S,hyp.a','lower',otherOptions{:});
            ub = supportFunc_(S,hyp.a','upper',otherOptions{:});

            if ~contains_(interval(lb,ub),hyp.b,'exact',0)
               res = false;
               return;
            end

            % check intersection with inequality constraints
            C = hyp.C; d = hyp.d;

            for i = 1:size(C,1)
               bound = supportFunc_(S,C(i,:)','lower',otherOptions{:});
               if bound > d(i)
                  res = false;
                  return;
               end
            end
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_intersectPolyPoly(hyp,P)
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
A = [P.A;hyp.C];
b = [P.b;hyp.d];

% construct matrices for equality constraint
Aeq = hyp.a;
beq = hyp.b;

% introduce slack variables y
m = length(b);
n = size(A,2);

problem.Aineq = [A,-eye(m);zeros(m,n),-eye(m)];
problem.bineq = [b;zeros(m,1)];

problem.Aeq = [Aeq,zeros(1,m)];
problem.beq = beq;

problem.f = [zeros(n,1);ones(m,1)];

problem.lb = [];
problem.ub = [];

% solve the dual problem using linear programming
[~,val,exitflag] = CORAlinprog(problem);

% check if intersection between the two polytopes is empty
res = true;

if exitflag < 0 || val > eps
    res = false; 
end

end
    
function res = aux_intersectPolyZonoBundle(obj1,obj2)
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
h = obj1.a';
d = obj1.b;

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
   
    Z = obj2.Z{i};
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

temp = repmat(eye(n),[obj2.parallelSets,1]);
Aeq = [temp,zeros(size(temp,1),p),Aeq];

% construct objective function
f = zeros(size(Aeq,2),1);
f(n+1:n+p) = ones(p,1);

problem.f = f;
problem.Aineq = A;
problem.bineq = b;
problem.Aeq = Aeq;
problem.beq = beq;
problem.lb = [];
problem.ub = [];

% solve linear program
[~,val,exitflag] = CORAlinprog(problem);

% check if intersection between the two polytopes is empty
res = true;

tol = 1e-8;
if exitflag < 0 || (val > 0 && ~withinTol(val,0,tol))
    res = false; 
end

end

% ------------------------------ END OF CODE ------------------------------
