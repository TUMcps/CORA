function [res,msg] = c_V(val,sys,params,options)
% c_V - costum validation function to check whether
%    params.V is of right size (different depending on system class)
%
% Syntax:
%    [res,msg] = c_V(val,sys,params,options)
%
% Inputs:
%    val - value for options.nrConstInp
%    sys - some contDynamics object
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
% Written:       10-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

if isa(sys,'linearSys') || isa(sys,'linearSysDT')
    if dim(val) ~= sys.nrOfOutputs
        % has to be equal to the number of system outputs
        res = false;
        msg = 'has to be equal to the number of outputs';
    end
else
    % nonlinearSysDT: currently no check... integrate once output equations
    % are integrated in nonlinearSys classes
end

end

% ------------------------------ END OF CODE ------------------------------
