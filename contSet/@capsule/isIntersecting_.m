function res = isIntersecting_(C,S,type,tol,varargin)
% isIntersecting_ - determines if a capsule intersects a set
%
% Syntax:
%    res = isIntersecting_(C,S,type,tol)
%
% Inputs:
%    C - capsule object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    C1 = capsule([0;0],[-2;2],2);
%    C2 = capsule([2;2],[-0.5;1],1);
%    C3 = capsule([3;3],[-0.5;1],1);
%
%    isIntersecting(C1,C2)
%    isIntersecting(C1,C3)
%
%    figure; hold on
%    plot(C1,[1,2],'b');
%    plot(C2,[1,2],'g');
%
%    figure; hold on
%    plot(C1,[1,2],'b');
%    plot(C3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, conZonotope/isIntersecting_

% Authors:       Niklas Kochdumper
% Written:       21-November-2019 
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[C,S] = reorderNumeric(C,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < C.precedence
    res = isIntersecting_(S,C,type,tol);
    return
end

% numeric case: check containment
if isnumeric(S)
    res = contains_(C,S,type,tol);
    return
end

% sets must not be empty
if representsa_(C,'emptySet',0) || representsa_(S,'emptySet',0)
    res = false;
    return
end

% capsule-capsule intersection check via shortest distance
if isa(S,'capsule')
    res = aux_isIntersecting_capsule(C,S);
    return
end

if isa(S,'polytope') && representsa_(S,'hyperplane',1e-12)
    res = aux_isIntersecting_hyperplane(C,S,type,tol);
    return
end

% for other sets, no exact algorithm supported
if isa(S,'contSet') && strcmp(type,'exact')
    throw(CORAerror('CORA:noExactAlg',C,S));
end

% approximate algorithm
if isa(S,'contSet') && strcmp(type,'approx')

    % enclose by zonotope and try again
    Z_C = zonotope(C);
    if S.precedence > Z_C.precedence
        res = isIntersecting_(S,Z_C,type,tol); 
    else
        res = isIntersecting_(Z_C,S,type,tol); 
    end
    return
end

throw(CORAerror('CORA:noops',C,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isIntersecting_capsule(C1,C2)
% check if two capsules C1 and C2 intersect

    % get object properties
    c1 = C1.c;
    c2 = C2.c;
    
    n = length(c1);
    
    g1 = C1.g;
    if isempty(g1)
        g1 = zeros(n,1);
    end
    
    g2 = C2.g;
    if isempty(g2)
        g2 = zeros(n,1); 
    end
    
    % construct objective function and constraints for quadratic program
    problem.H = 0.5*[g1'*g1 -g1'*g2; -g1'*g2 g2'*g2];
    problem.f = 2*[(c1-c2)'*g1;-(c1-c2)'*g2];
    problem.Aineq = [];
    problem.bineq = [];
    problem.Aeq = [];
    problem.beq = [];
    problem.lb = [-1;-1];
    problem.ub = [1;1];
    
    % solve quadratic program
    x = CORAquadprog(problem);
    
    % compute minimum distance between point on the axis
    p1 = c1 + g1*x(1);
    p2 = c2 + g2*x(2);
    
    dist = norm(p1-p2);
    
    % check for intersection
    tmp = C1.r + C2.r;
    res = dist < tmp | withinTol(dist,tmp);
    
end

function res = aux_isIntersecting_hyperplane(C,P,type,tol)

res = true;

% check intersection with hyperplane
I = supportFunc_(C,P.Ae(1,:)','range');
if ~contains_(I,P.be(1),'exact',tol)
    res = false;
    return
end

end

% ------------------------------ END OF CODE ------------------------------
