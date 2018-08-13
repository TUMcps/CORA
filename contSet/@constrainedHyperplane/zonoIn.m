function res = zonoIn(obj, Z)
% zonoIn - Checks if a zonotope is within a halfspace
%
% Syntax:  
%    res = zonoIn(obj, Z)
%
% Inputs:
%    obj - halfspace object
%    Z - zonotope object
%
% Outputs:
%    res - boolean variable, 1 if in halfspace, 0 otherwise
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      10-March-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute if hyperplane intersects
planeInter = zonoIn(obj.h, Z);

%check if constraints are satisfied
cSat = constrSat(Z, obj.C, obj.d);

%return overall result
res = planeInter & cSat;


%------------- END OF CODE --------------
