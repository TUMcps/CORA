function res = mid(obj)
% mid - returns the center of an interval
%
% Syntax:  
%    res = mid(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    res - numerical value
%
% Example: 
%    a = interval([-1 1], [1 2]);
%    b = mid(a)
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

res = 0.5*(obj.inf + obj.sup);

%------------- END OF CODE --------------