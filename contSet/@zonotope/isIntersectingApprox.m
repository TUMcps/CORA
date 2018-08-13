function res = isIntersectingApprox(obj1,obj2)
% isIntersectingApprox - approximate test if the zonotope obj1 intersects 
%                        obj2. If a intersection occurs, the function
%                        always returns 1. But the function possibly also
%                        returns 1 if no intersection occurs
%
% Syntax:  
%    res = isIntersectingApprox(obj1,obj2)
%
% Inputs:
%    obj1 - zonotope object
%    obj2 - interval or zonotope
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
% Written:      24-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% different set representations
if isa(obj2,'zonotope')
    
    % convert zonotopes to intervals and check for intersection
    res = isIntersecting(interval(obj1),interval(obj2));
    
elseif isa(obj2,'interval')
    
    % convert zonotopes to intervals and check for intersection
    res = isIntersecting(obj1,interval(obj2));
    
else
    error('Operation not implemented yet!');
end


%------------- END OF CODE --------------