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
% Written:      02-September-2013
% Last update:  22-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%compute map on normal vector of the halfspace
IH = interval(obj.c' * Z) + (-obj.d);

%check if interval contains 0
if supremum(IH)<=0
    res = 1;
else
    res = 0;
end


%------------- END OF CODE --------------