function res = isIntersectingApprox(obj1,obj2)
% isIntersectingApprox - approximate test if a constrained hyperplane obj1 intersects 
%                        obj2. If a intersection occurs, the function
%                        always returns 1. But the function possibly also
%                        returns 1 if no intersection occurs
%
% Syntax:  
%    res = isIntersectingApprox(obj1,obj2)
%
% Inputs:
%    obj1 - constrainedHyperplane object
%    obj2 - interval, zonotope or zonotopeBundle object
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
% Written:      16-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


if isa(obj2,'zonotopeBundle')
   
    res = 1;
    
    % if one of the parallel sets does not intersect the hyperplane, then
    % the zonotopeBundle (=intersection of sets) does not intersect the 
    % hyperplane either
    for i = 1:length(obj2.Z)
       if ~isIntersectingApprox(obj1,obj2.Z{i})
          res = 0;
          return;
       end
    end
    
    
elseif isa(obj2,'zonotope') || isa(obj2,'interval')
   
    % check for intersection with hyperplane
    res = isIntersecting(obj1.h,obj2);
    
    % check if the set is located fully outside the feasible region defined
    % by the constraints
    if res
        
       % extract set midpoint
       if isa(obj2,'zonotope')
           c = center(obj2);
       else
           c = mid(obj2); 
       end      
        
       % loop over all constraints
       for i = 1:length(obj1.d)
           
           % check if center outside feasible region
           if obj1.C(i,:)*c > obj1.d(i)
               
               % center outside + set does not intersect constraint 
               % hyperplane -> set fully outside feasible region
               temp = halfspace(obj1.C(i,:),obj1.d(i));
               
               if ~isIntersecting(temp,obj2)
                   res = 0;
                   return;
               end
           end        
       end
    end
    
else
    error('Function not implemented for this kind of sets!');
end


%------------- END OF CODE --------------