function res = rad(obj)
% rad - returns the radius of an interval
%
% Syntax:  
%    res = rad(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    res - numerical value
%
% Example: 
%    a = interval([-1 1], [1 2]);
%    b = rad(a)
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

res = 0.5*(obj.sup - obj.inf);

%------------- END OF CODE --------------