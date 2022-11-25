function res = infimum(obj)
% infimum - returns the infimum of an interval
%
% Syntax:  
%    res = inf(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    res - numerical value
%
% Example: 
%    a = interval([-1 1], [1 2]);
%    b = inf(a)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      25-June-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = obj.inf;

%------------- END OF CODE --------------