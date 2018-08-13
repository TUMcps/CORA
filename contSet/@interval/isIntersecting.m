function res = isIntersecting(obj1,obj2)
% isIntersecting - determines if a set obj2 intersects an interval obj1
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%
% Inputs:
%    obj1 - interval object
%    obj2 - some set
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

% Author:       Matthias Althoff
% Written:      22-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%make an intervalhull of obj2
I = interval(obj2);

%check if intervals intersect
int = (supremum(obj1)<infimum(I))|(supremum(I)<infimum(obj1));
res = ~any(int);

%------------- END OF CODE --------------