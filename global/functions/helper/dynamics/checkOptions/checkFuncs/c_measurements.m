function [res,msg] = c_measurements(val,sys,params)
% c_measurements - costum validation function to check whether params.y
%    and obj.dt match
%
% Syntax:
%    [res,msg] = c_measurements(val,sys,params,options)
%
% Inputs:
%    val - value for params.y
%    sys - some contDynamics object
%    params - model parameters
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
% Written:       08-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

if isa(sys,'linearSysDT') || isa(sys,'nonlinearSysDT')
    % check if size of params.y matches options.timeStep
    if (size(val,2) == 1 && ~any(val))
        % params.y is default value (zero-vector)
        return;
    elseif size(val,2) ~= round((params.tFinal-params.tStart) / sys.dt)
        res = false;
        msg = ['does not comply with obj.dt: \n'...
            'The number of steps in the reachability analysis given by \n'...
            '   (params.tFinal-params.tStart)/obj.dt\n'...
            'has to match the number of columns in params.y'];
        return;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
