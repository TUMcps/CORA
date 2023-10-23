function res = isIntersecting_(C,S,varargin)
% isIntersecting_ - determines if a capsule intersects a set
%
% Syntax:
%    res = isIntersecting_(C,S)
%    res = isIntersecting_(C,S,type)
%
% Inputs:
%    C - capsule object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
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
    
    % check for intersection
    if isa(S,'capsule')
        res = aux_intersectionCapsule(C,S);

    elseif isa(S,'conHyperplane')
        res = isIntersecting_(S,C,type);

    else
        % exact or over-approximative algorithm
        if strcmp(type,'exact')
            throw(CORAerror('CORA:noExactAlg',C,S));
        else
            res = isIntersecting_(zonotope(C),S,type); 
        end
    end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_intersectionCapsule(C1,C2)
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
    H = 0.5*[g1'*g1 -g1'*g2; -g1'*g2 g2'*g2];
    f = 2*[(c1-c2)'*g1;-(c1-c2)'*g2];
    lb = [-1;-1];
    ub = [1;1];
    
    % solve quadratic program
    options = optimoptions('quadprog','display','off');
    
    x = quadprog(H,f,[],[],[],[],lb,ub,[0;0],options);
    
    % compute minimum distance between point on the axis
    p1 = c1 + g1*x(1);
    p2 = c2 + g2*x(2);
    
    dist = norm(p1-p2);
    
    % check for intersection
    tmp = C1.r + C2.r;
    res = dist < tmp | withinTol(dist,tmp);
    
end

% ------------------------------ END OF CODE ------------------------------
