function res = acos( obj )
% acos - Overloaded 'acos()' operator for a taylm expression
%
% Syntax:  
%    res = acos(obj)
%
% Inputs:
%    obj - a taylm expression
%
% Outputs:
%    res - a taylm expression object
%
% Example: 
%
% Other m-files required: taylm, interval
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, asin
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Dmitry Grebenyuk
% Written:      19-August-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------
    res = pi/2 - asin(obj);
end
%------------- END OF CODE --------------