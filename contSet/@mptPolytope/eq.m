function res = eq(obj, obj2)
% eq - overloads the '==' operator,
% returns 1 if two polytopes are equal and 0 otherwise
%
% Syntax:  
%    res = eq(obj, obj2)
%
% Inputs:
%    obj - mptPolytope object
%    obj2 - mptPolytope object
%
% Outputs:
%    res - boolean whether equal or not
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
% Written:      31-August-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = (obj.P == obj2.P);

%------------- END OF CODE --------------