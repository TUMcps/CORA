function res = isIntersecting(obj1,obj2,varargin)
% isIntersecting - determines if capsule obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%    res = isIntersecting(obj1,obj2,type)
%
% Inputs:
%    obj1 - capsule object
%    obj2 - conSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    C1 = capsule([0;0],[-2;2],2);
%    C2 = capsule([2;2],[-0.5;1],1);
%    C3 = capsule([3;3],[-0.5;1],1);
%
%    isIntersecting(C1,C2)
%    isIntersecting(C1,C3)
%
%    figure
%    hold on
%    plot(C1,[1,2],'b');
%    plot(C2,[1,2],'g');
%
%    figure
%    hold on
%    plot(C1,[1,2],'b');
%    plot(C3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isIntersecting

% Author:       Niklas Kochdumper
% Written:      21-Nov-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && ~isempty(varargin{1}) 
        type = varargin{1};
    end

    % get capsule object
    if ~isa(obj1,'capsule')
        temp = obj1;
        obj1 = obj2;
        obj2 = temp;
    end
    
    % check for intersection
    if isa(obj2,'capsule')
        % check for dimension mismatch
        if dim(obj1) ~= dim(obj2)
            [id,msg] = errDimMismatch();
            error(id,msg);
        end
        res = intersectionCapsule(obj1,obj2);
    elseif isa(obj2,'hyperplane')
        res = isIntersecting(obj2,obj1,type); 
    else
       
       % exact or over-approximative algorithm
       if strcmp(type,'exact')
           error('No exact algorithm implemented for this set representation!');
       else
           res = isIntersecting(zonotope(obj1),obj2); 
       end     
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = intersectionCapsule(C1,C2)
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
    res = dist <= C1.r + C2.r;
    
    % additional check for closeness of distances
    if ~res
        dummy = C1.r + C2.r; dummy(dummy==0) = 1;
        if abs(dist-dummy)/dummy < 1e-9
            % tight enough -> intersection at one point
            res = true;
        end
    end
    
end

%------------- END OF CODE --------------