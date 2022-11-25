function res = le(obj, obj2)
% le - returns 1 if one polytopes is equal or enclosed by the other one and 
% 0 otherwise; overloads the <= operator
%
% Syntax:  
%    res = le(obj, obj2)
%
% Inputs:
%    obj - mptPolytope object
%    obj2 - mptPolytope object
%
% Outputs:
%   res - result in {0,1}
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
% Written:      14-November-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = (obj.P <= obj2.P);

%------------- END OF CODE --------------