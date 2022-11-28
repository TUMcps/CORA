function res = contains(C,S)
% contains - determines if a capsule contains a set or a point
%
% Syntax:  
%    res = contains(C,S)
%
% Inputs:
%    C - capsule object
%    S - contSet object or single point
%
% Outputs:
%    res - true/false
%
% Example: 
%    C1 = capsule([0;0],[-2;2],2);
%    C2 = capsule([-0.5;0],[-0.5;1],1);
%    C3 = capsule([-2;0],[-0.5;1],1);
%
%    contains(C1,C2)
%    contains(C1,C3)
%
%    figure; hold on;
%    plot(C1,[1,2],'b');
%    plot(C2,[1,2],'g');
%
%    figure; hold on;
%    plot(C1,[1,2],'b');
%    plot(C3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/contains

% Author:       Niklas Kochdumper
% Written:      20-November-2019
% Last update:  15-November-2022 (MW, return logical array for points)
%               25-November-2022 (MW, rename 'contains')
% Last revision:---

%------------- BEGIN CODE --------------
    
    % pre-processing
    [res,vars] = pre_contains('capsule',C,S);
    
    % check premature exit
    if res
        % if result has been found, it is stored in the first entry of var
        res = vars{1}; return
    end

    
    % init result
    res = true;

    % point in capsule containment
    if isnumeric(S)
        
        res = false(1,size(S,2));
        for i = 1:size(S,2)
            res(i) = containsPoint(C,S(:,i));
        end
        
    % capsule in capsule containment    
    elseif isa(S,'capsule')
        
        res = containsCapsule(C,S);
        
    % non polytopic set in capsule containment
    elseif isa(S,'ellipsoid') || isa(S,'taylm') || isa(S,'polyZonotope')
       
        res = contains(C,zonotope(S));
        
    % polytopic set in capsule containment    
    else
        
        % compute vertices
        V = vertices(S);
        
        % check containment for each vertex
        for i = 1:size(V,2)
            if ~containsPoint(C,V(:,i))
                res = false;
                return;
            end
        end     
    end
    
end


% Auxiliary Functions -----------------------------------------------------

function res = containsCapsule(C1,C2)
% checks if the capsule obj2 is contained in the capsule obj2

    if isempty(C2.g)
       res = containsSphere(C1,C2.c,C2.r); 
    else
       res1 = containsSphere(C1,C2.c+C2.g,C2.r);
       res2 = containsSphere(C1,C2.c-C2.g,C2.r);
       
       res = res1 & res2;
    end
end

function res = containsSphere(C,c,r)
% checks if the capsule contains the hypersphere defined by the center 
% and the radius
    
    % check case where capsule is a hypersphere (no generator)
    if isempty(C.g)
        
        tmp = norm(C.c-c) + r;
        res = tmp < C.r | withinTol(tmp,C.r);
       
    else
       
        ng = norm(C.g);
        g_ = C.g/ng;
        proj = (c-C.c)'*g_;
        
        % check if center is in hypercylinder
        if abs(proj) < norm(C.g)
           
            % compute distance to axis
            diff = (C.c + proj*g_) - c;
            
            % check if distance to axis is smaller than the radius
            tmp = norm(diff) + r;
            res = tmp < C.r | withinTol(tmp,C.r);
            
        else    
            % check if point is in upper or lower hypersphere
            tmp = norm(C.c + C.g - c) + r;
            res1 = tmp < C.r | withinTol(tmp,C.r);
            tmp = norm(C.c - C.g - c) + r;
            res2 = tmp < C.r | withinTol(tmp,C.r);
            
            res = res1 | res2;
        end  
    end
end

function res = containsPoint(C,p)
% checks if a point is contained in the capsule
   
    % get object properties
    g = C.g;
    r = C.r;
    c = C.c;

    % check case where capsule is a hypersphere
    if isempty(g)
        
       res = inSphere(c,r,p);
       
    else
       
        ng = norm(g);
        g_ = g/ng;
        proj = (p-c)'*g_;
        
        % check if point is in hypercylinder
        if abs(proj) < norm(g)
           
            % compute distance to axis
            diff = (c + proj*g_) - p;
            
            % check if distance to axis is smaller than the radius
            tmp = norm(diff);
            res = tmp < r | withinTol(tmp,r);
            
        else    
            % check if point is in upper or lower hypersphere
            res = inSphere(c+g,r,p) | inSphere(c-g,r,p);
        end  
    end
end


function res = inSphere(c,r,p)
% checks if a point is contained in the hypersphere defined by center and
% radius
    
    tmp = norm(p-c);
    res = tmp < r | withinTol(tmp,r);

end

%------------- END OF CODE --------------