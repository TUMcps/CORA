function obj = hull(intVal1, intVal2)
% hull - The union of two intervals
% a = hull(b, c);
%
% Syntax:  
%    obj = hull(intVal1, intVal2)
%
% Inputs:
%    intVal1, intVal2 - intervals 
%
% Outputs:
%    obj - an interval object 
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Dmitry Grebenyuk
% Written:      20-April-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

obj = intVal1;

obj.inf = min(intVal1.inf, intVal2.inf);
obj.sup = max(intVal1.sup, intVal2.sup);


%------------- END OF CODE --------------
