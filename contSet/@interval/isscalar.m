function res = isscalar(intVal)
% isscalar - returns 1 if the interval is scalar and 0 otherwise
%
% Syntax:  
%     res = isscalar(intVal)
%
% Inputs:
%    intVal - interval object
%
% Outputs:
%    res - bool
%
% Example: 
%    a = interval(-1,2);
%    res = isscalar(a)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff
% Written:      25-June-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain result; only infimum is checked for efficiency reasons
res = isscalar(intVal.inf);

%------------- END OF CODE --------------