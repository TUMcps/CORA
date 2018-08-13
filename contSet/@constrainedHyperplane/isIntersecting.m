function res = isIntersecting(obj1,obj2)
% isIntersecting - determines if constrained hyperplane obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%
% Inputs:
%    obj1 - constrainedHyperplane object
%    obj2 - interval, zonotope or constrained zonotope object
%
% Outputs:
%    result - 1/0 if set is intersecting, or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      22-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = 0;

% convert interval to zontope
if isa(obj2,'interval')
   obj2 = zonotope(obj2); 
end

% convert zonotope to constrained zonotope
if isa(obj2,'zonotope')
   obj2 = conZonotope(obj2.Z,[],[]); 
end

% different set representations
if isa(obj2,'conZonotope')
    
    % check if the unconstrained zonotope intersects the hyperplane
    zono = zonotope(obj2.Z);
    
    if isIntersecting(obj1.h,zono)

        % calculate intersection between hyperplane and constrained zonotope
        interZono = obj2 & obj1.h;

        % check if part of the intersection is located inside the feasible
        % region defined by the constraints
        for i = 1:size(obj1.C,1)
           bound = boundDir(interZono,obj1.C(i,:)','lower');
           if bound < obj1.d(i)
                res = 1;
                break;
           end
        end
    end
    
else
   
    error('Operation not implemented yet!');
    
end


%------------- END OF CODE --------------