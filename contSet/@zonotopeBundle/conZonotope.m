function res = conZonotope(obj)
% conZonotope - convert a zonotopeBundle object into a constrained 
%               zonotope object
%
% Syntax:  
%    res = conZonotope(obj)
%
% Inputs:
%    obj - zonotopeBundle object
%
% Outputs:
%    res - c-zonotope object
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      23-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
% initialization
res = conZonotope(obj.Z{1});

% calculate the intersection of the parallel sets
for i = 2:obj.parallelSets
   temp = conZonotope(obj.Z{i});
   res = res & temp;
end

%------------- END OF CODE --------------