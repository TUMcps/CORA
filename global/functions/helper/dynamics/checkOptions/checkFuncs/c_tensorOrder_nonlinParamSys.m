function [res,msg] = c_tensorOrder_nonlinParamSys(val,sys,params,options)
% c_tensorOrder_nonlinParamSys - costum validation function for
%    options.tensorOrder (only for class nonlinParamSys)
%
% Syntax:
%    [res,msg] = c_tensorOrder_nonlinParamSys(val,sys,params,options)
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
    if val ~= 3
        msg = 'has to be 3 if options.alg = ''poly''';
        res = false; return;
    end
    
% options.alg = 'lin' or 'linRem'    
elseif any(strcmp(options.alg,{'lin','linRem'}))

    if isa(params.paramInt,'interval')
        if val ~= 2
            msg = 'has to be 3 if params.paramInt is an interval';
            res = false; return;
        end 
    elseif isvector(params.paramInt)
        if ~any(val == [2,3])
            msg = 'has to be either 2 or 3';
            res = false; return;
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
