function [res,msg] = c_inputTraj(val,sys,params,options)
% c_inputTraj - costum validation function to check whether params.u
%    and options.timeStep match
%
% Syntax:
%    [res,msg] = c_inputTraj(val,sys,params,options)
%
% Inputs:
%    val - value for options.timeStep
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

% Author:       Mark Wetzlinger
% Written:      25-May-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% assume check ok
res = true;
msg = '';

% check if size of params.u matches options.timeStep
if ~isfield(params,'u') || (size(params.u,2) == 1 && ~any(params.u))
    % params.u not provided or default value (zero-vector)
    return;
elseif size(params.u,2) ~= round((params.tFinal-params.tStart) / options.timeStep)
    res = false;
    msg = ['does not comply with params.u: \n'...
        'The number of steps in the reachability analysis given by \n'...
        '   (params.tFinal-params.tStart)/options.timeStep\n'...
        'has to match the number of columns in params.u'];
    return;
end

end

%------------- END OF CODE --------------
