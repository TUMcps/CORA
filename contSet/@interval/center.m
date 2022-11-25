function res = center(obj)
% center - returns the center of an interval
%
% Syntax:  
%    res = center(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    res - center of interval (vector)
%
% Example: 
%    a = interval([-1 1], [1 2]);
%    b = center(a)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      26-June-2015
% Last update:  02-Sep-2019 (rename mid -> center)
% Last revision:---

%------------- BEGIN CODE --------------

res = 0.5*(obj.inf + obj.sup);

%------------- END OF CODE --------------