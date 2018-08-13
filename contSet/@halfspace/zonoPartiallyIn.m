function res = zonoPartiallyIn(obj, Z)
% zonoIn - Checks if a zonotope is partially within a halfspace
%
% Syntax:  
%    res = zonoPartiallyIn(obj, Z)
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
% Written:      02-September-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute map on normal vector of the halfspace
IH = interval(obj.c' * Z) + (-obj.d);

%check if interval contains 0
if IH(1,1)<=0
    res = 1;
else
    res = 0;
end


%------------- END OF CODE --------------