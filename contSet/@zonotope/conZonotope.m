function res = conZonotope(obj)
% conZonotope - convert a zonotope object into a constrained zonotope
%               object
%
% Syntax:  
%    res = conZonotope(obj)
%
% Inputs:
%    obj - zonotope object
%
% Outputs:
%    res - c-zonotope object
%
% Example:
%    Z = zonotope([1;0],[1 0; 0 1]);
%    cZ = conZonotope(Z);
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
    
res = conZonotope(obj.Z);

%------------- END OF CODE --------------