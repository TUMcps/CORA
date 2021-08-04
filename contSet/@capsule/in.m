function res = in(obj1,obj2)
% in - determines if obj2 is fully inside a halfspace obj1
%
% Syntax:  
%    res = in(obj1,obj2)
%
% Inputs:
%    obj1 - capsule object
%    obj2 - contSet object or single point
%
% Outputs:
%    res - 1/0 if obj2 is contained in obj1, or not
%
% Example: 
%    C1 = capsule([0;0],[-2;2],2);
%    C2 = capsule([-0.5;0],[-0.5;1],1);
%    C3 = capsule([-2;0],[-0.5;1],1);
%
%    in(C1,C2)
%    in(C1,C3)
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
% See also: zonotope/in

% Author:       Niklas Kochdumper
% Written:      20-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = 1;

    % point in capsule containment
    if isnumeric(obj2)
        
        for i = 1:size(obj2,2)
            res = containsPoint(obj1,obj2(:,i));
            if res ~= 1
                return;
            end
        end
        
    % capsule in capsule containment    
    elseif isa(obj2,'capsule')
        
        % check dimension mismatch
        if dim(obj1) ~= dim(obj2)
            [id,msg] = errDimMismatch();
            error(id,msg);
        end
        
        res = containsCapsule(obj1,obj2);
        
    % non polytopic set in capsule containment
    elseif isa(obj2,'ellipsoid') || isa(obj2,'taylm') || isa(obj2,'polyZonotope')
       
        res = in(obj1,zonotope(obj2));
        
    % polytopic set in capsule containment    
    else
        
        % compute vertices
        V = vertices(obj2);
        
        % check containment for each vertex
        for i = 1:size(V,2)
            if ~containsPoint(obj1,V(:,i))
                res = 0;
                return;
            end
        end     
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = containsCapsule(obj1,obj2)
% checks if the capsule obj2 is contained in the capsule obj2

    if isempty(obj2.g)
       res = containsSphere(obj1,obj2.c,obj2.r); 
    else
       res1 = containsSphere(obj1,obj2.c+obj2.g,obj2.r);
       res2 = containsSphere(obj1,obj2.c-obj2.g,obj2.r);
       
       res = res1 & res2;
    end
end

function res = containsSphere(obj,c,r)
% checks if the capsule obj contains the hyperspere defined by the center c
% and the radius r
    
    % check case where capsule is a hyobj.c + proj*obj.gpersphere
    if isempty(obj.g)
        
       res = norm(obj.c-c) + r <= obj.r;
       
    else
       
        ng = norm(obj.g);
        g_ = obj.g/ng;
        proj = (c-obj.c)'*g_;
        
        % check if center is in hypercylinder
        if abs(proj) < norm(obj.g)
           
            % compute distance to axis
            diff = (obj.c + proj*g_) - c;
            
            % check if distance to axis is smaller than the radius
            res = norm(diff) + r <= obj.r;
            
        else    
            % check if point is in upper or lower hypersphere
            res1 = norm(obj.c+obj.g-c) + r <= obj.r;
            res2 = norm(obj.c-obj.g-c) + r <= obj.r;
            
            res = res1 | res2;
        end  
    end
end

function res = containsPoint(obj,p)
% checks if the point p is contained by the capsule obj
   
    % get object properties
    g = obj.g;
    r = obj.r;
    c = obj.c;

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
            res = norm(diff) <= r;
            
        else    
            % check if point is in upper or lower hypersphere
            res = inSphere(c+g,r,p) | inSphere(c-g,r,p);
        end  
    end
end


function res = inSphere(c,r,p)
% checks if the point p is in the hypersphere defined by center c and
% radius r
    
    res = norm(p-c) <= r;

end

%------------- END OF CODE --------------