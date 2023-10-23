function [res,msg] = c_scaleFac(val,sys,options)
% c_scaleFac - costum validation function for options.factor
%
% Syntax:
%    [res,msg] = c_scaleFac(val,sys,list)
%
% Inputs:
%    val - value for given param / option
%    sys - linearSys object
%    options - algorithm parameters
%
% Outputs:
%    res - logical whether validation was successful
%    msg - error message if validation failed
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Authors:       Mark Wetzlinger
% Written:       08-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

if ischar(val)
    if ~strcmp(val,'auto')
        res = false;
    end
elseif isnumeric(val) && isscalar(val)
	if val <= 0 || val > 1
        res = false;
    end
else
    res = false;
end

% unified error message for all cases
if ~res
    msg = 'has to be either ''auto'' or a scalar in [0,1]';
end


end

% ------------------------------ END OF CODE ------------------------------
