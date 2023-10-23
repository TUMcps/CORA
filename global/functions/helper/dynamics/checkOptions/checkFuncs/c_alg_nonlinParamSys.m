function [res,msg] = c_alg_nonlinParamSys(val,sys,params,options)
% c_alg_nonlinParamSys - costum validation function for
%    options.tensorOrder (only for class nonlinParamSys)
%
% Syntax:
%    [res,msg] = c_alg_nonlinParamSys(val,sys,params,options)
%
% Inputs:
%    val - value for given param / option
%    sys - contDynamics object
%    params - model parameters
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
% Written:       03-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

% options.alg = 'poly'
if strcmp(options.alg,'poly')
    if ~isvector(params.paramInt)
        msg = 'can only be ''poly'' if params.paramInt is a vector';
        res = false; return;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
