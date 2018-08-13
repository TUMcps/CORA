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
% Written:      08-August-2011
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%compute map on normal vector of the halfspace
IH = interval(obj.c' * Z) + (-obj.d);

%check if interval contains 0
if infimum(IH)<0 && supremum(IH)>=0
    res = 1;
else
    res = 0;
end


%------------- END OF CODE --------------