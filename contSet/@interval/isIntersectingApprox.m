function res = isIntersectingApprox(obj1,obj2)
% isIntersectingApprox - approximate test if the interval obj1 intersects 
%                        obj2. If a intersection occurs, the function
%                        always returns 1. But the function possibly also
%                        returns 1 if no intersection occurs
%
% Syntax:  
%    res = isIntersectingApprox(obj1,obj2)
%
% Inputs:
%    obj1 - interval object
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
    
elseif isa(obj2,'interval')
   
    % caclulate exact intersection
    res = isintersecting(obj1,obj2);
    
elseif isa(obj2,'zonotope')
   
    % caclulate exact intersection
    res = isIntersecting(obj1,interval(obj2));
    
else
    error('Function not implemented for this kind of sets!');
end


%------------- END OF CODE --------------