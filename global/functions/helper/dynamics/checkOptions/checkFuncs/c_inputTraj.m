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

% Authors:       Mark Wetzlinger
% Written:       25-May-2021
% Last update:   17-November-2021 (MW, allowed different lengths depending on existence of feedthrough matrix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

% check if size of params.u matches options.timeStep
if ~isfield(params,'u') || (size(params.u,2) == 1 && ~any(params.u))
    % params.u not provided or default value (zero-vector)
    return;
else
    steps = round((params.tFinal-params.tStart) / val);
    if (isa(sys,'linearSys') || isa(sys,'linearSysDT')) && any(any(sys.D))
        % feedthrough matrix given, for last output we require one more
        % entry in params.u
        if size(params.u,2) ~= steps + 1
            res = false;
            msg = ['does not comply with params.u:\n'...
                'Case: Feedthrough matrix D provided\n'...
                'The number of steps in the reachability analysis plus 1 given by\n'...
                '   (params.tFinal-params.tStart)/options.timeStep + 1\n'...
                'has to match the number of columns in params.u'];
            return;
        end
    else
        % feedthrough matrix is all-zero
        if ~any(size(params.u,2) == [steps, steps+1])
            res = false;
            msg = ['does not comply with params.u:\n'...
                'Case: Not a linear system / no feedthrough matrix provided\n'...
                'The number of steps in the reachability analysis given by\n'...
                '   (params.tFinal-params.tStart)/options.timeStep\n'...
                'has to match the number of columns in params.u'];
            return;
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
