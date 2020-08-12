function res = conZonotope(obj)
% conZonotope - Converts an interval object into a conZonotope object
%
% Syntax:  
%    res = conZonotope(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    res - conZonotope object
%
% Example: 
%    I = interval([1;-1], [2; 1]);
%    cZ = conZonotope(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Author:       Niklas Kochdumper
% Written:      21-Nov-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = conZonotope(zonotope(obj));

%------------- END OF CODE --------------