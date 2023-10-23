function res = getCoef(obj)
% getCoef - returns coefficients
%
% Syntax:
%    res = getCoef(obj)
%
% Inputs:
%    obj - a Taylor model
%
% Outputs:
%    res - array of numbers 
%
% Example: 
%
% Other m-files required: taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Dmitry Grebenyuk
% Written:       06-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%	res = arrayfun(@(a) s_getCoef(a), obj, 'UniformOutput', 0);
%
%end
%
%function res = s_getCoef( obj )
%
res = obj.coefficients;

% ------------------------------ END OF CODE ------------------------------
