function [res,msg] = c_tensorOrder(val,sys,options)
% c_tensorOrder - costum validation function for options.tensorOrder
%
% Syntax:
%    [res,msg] = c_tensorOrder(val,sys,options)
%
% Inputs:
%    val - value for given param / option
%    sys - contDynamics object
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
% Written:       04-March-2019
% Last update:   21-April-2020 (split in lin/poly)
%                03-May-2020 (rewriting of error msgs using class(obj))
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

% set min / max admissible values depending on other information
maxOrder = 7;
if isa(sys,'nonlinearSysDT')
    minOrder = 2;
elseif strcmp(options.alg,'lin')
    minOrder = 2;
elseif strcmp(options.alg,'linRem')
    minOrder = 2;
    maxOrder = 2;
elseif strcmp(options.alg,'poly')
    minOrder = 3;
    if isfield(options,'approxDepOnly') && options.approxDepOnly
        minOrder = 2;
    end
end
% perform check
if val < minOrder || val > maxOrder
    res = false;
    msg = ['has to be between ' num2str(minOrder) ' and ' num2str(maxOrder) ...
        ' for chosen options.alg = ' options.alg];
end

end

% ------------------------------ END OF CODE ------------------------------
