function [res,msg] = c_inputTrajDT(val,sys,params)
% c_inputTrajDT - costum validation function to check whether params.u
%    and obj.dt match
%
% Syntax:
%    [res,msg] = c_inputTrajDT(val,sys,params,options)
%
% Inputs:
%    val - value for params.u
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
% Written:       25-May-2021
% Last update:   17-November-2021 (MW, different lengths allowed depending on existence of feedthrough matrix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

% check if size of params.u matches sys.dt
if (size(val,2) == 1 && ~any(val))
    % params.u is default value (zero-vector)
    return;
else
    steps = round((params.tFinal-params.tStart) / sys.dt);
    if isa(sys,'linearSysDT') && any(any(sys.D))
        % feedthrough matrix given, for last output we require one more
        % entry in params.u
        if size(val,2) ~= steps + 1
            res = false;
            msg = ['does not comply with obj.dt:\n'...
                'Case: Feedthrough matrix D provided\n'...
                'The number of steps in the reachability analysis plus 1 given by\n'...
                '   (params.tFinal-params.tStart)/obj.dt + 1\n'...
                'has to match the number of columns in params.u'];
            return;
        end
    else
        % feedthrough matrix is all-zero
        if ~any(size(val,2) == [steps, steps+1])
            res = false;
            msg = ['does not comply with obj.dt:\n'...
                'Case: No feedthrough matrix provided\n'...
                'The number of steps in the reachability analysis given by\n'...
                '   (params.tFinal-params.tStart)/obj.dt\n'...
                'has to match the number of columns in params.u'];
            return;
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
