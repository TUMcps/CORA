function res = zonoIntersect(obj, Z)
% zonoIntersect - Checks if a zonotope intersects with the halfspace
%
% Syntax:  
%    zonoIntersect(obj)
%
% Inputs:
%    obj - halfspace object
%
% Outputs:
%    ---
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      10-August-2011
% Last update:  14-May-2017
% Last revision:---

%------------- BEGIN CODE --------------

%compute if hyperplane intersects
planeInter = zonoIntersect(obj.h, Z);

%check if constraints are satisfied
cSat = constrSat_some(Z, obj.C, obj.d);

%return overall result
res = planeInter & cSat;


%------------- END OF CODE --------------