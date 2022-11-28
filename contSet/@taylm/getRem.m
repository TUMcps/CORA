function rem = getRem(obj)
% getRem - returns remainder
%
% Syntax:  
%    rem = getRem(obj)
%
% Inputs:
%    obj - a Taylor model
%
% Outputs:
%    rem - an interval 
%
% Example:
%    -
%
% Other m-files required: taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Dmitry Grebenyuk
% Written:      06-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
%
%	rem = arrayfun(@(a) s_getRem(a), obj, 'UniformOutput', 0);
%
%end
%
%function rem = s_getRem( obj )
%
    rem = obj.remainder;

%------------- END OF CODE --------------
