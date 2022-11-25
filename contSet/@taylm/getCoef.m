function res = getCoef( obj )
% getCoef( obj ) - returns coefficients
%
% Syntax:  
%    res = coef = getCoef( obj )
%
% Inputs:
%    obj - a Taylor model
%
% Outputs:
%    coef - array of numbers 
%
% Example: 
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

%	res = arrayfun(@(a) s_getCoef(a), obj, 'UniformOutput', 0);
%
%end
%
%function res = s_getCoef( obj )
%
    res = obj.coefficients;

end

