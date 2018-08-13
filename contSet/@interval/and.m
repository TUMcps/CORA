function obj = and(obj,otherObj)
% and - computes intersection of intervals. Overloaded '&' operator for
% intervals
%
% Syntax:  
%    obj = and(obj,otherObj)
%
% Inputs:
%    obj - interval object
%    otherObj - interval object
%
% Outputs:
%   obj - interval object
%
% Example: 
%    a = interval([1;-1], [2; 1]);
%    b = interval([1.5; -2], [2.5; 0]);
%    c = a&b
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      26-June-2015
% Last update:  --- 
% Last revision:---

%------------- BEGIN CODE --------------


inf = max(obj.inf, otherObj.inf);
sup = min(obj.sup, otherObj.sup);

if all(all(inf <= sup)) %twice all because of interval matrices
    obj.inf = inf;
    obj.sup = sup;
else
    obj.inf = [];
    obj.sup = [];
end



%------------- END OF CODE --------------