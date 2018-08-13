function res = isIntersecting(obj1,obj2)
% isIntersecting - determines if a set obj2 intersects a zonotope obj1
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%
% Inputs:
%    obj1 - zonotope object
%    obj2 - zonotope or interval object
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

% convert intervals to zonotopes
if isa(obj2,'interval')
   obj2 = zonotope(obj2); 
end

% different set representations
if isa(obj2,'zonotope')
    
    % convert to constraint zonotopes and check if the intersection is
    % empty
    res = ~isempty(conZonotope(obj1) & conZonotope(obj2));
    
else
    error('Operation not implemented yet!');
end

%------------- END OF CODE --------------