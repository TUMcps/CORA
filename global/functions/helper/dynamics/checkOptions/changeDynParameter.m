function options = changeDynParameter(errmsgid,params,options)
% changeDynParameter - Change dynamic parameter according to errmsg
%
% Syntax:
%    options = changeDynParameter(errmsgid,params,options)
%
% Inputs:
%    errmsgid - error message identifier
%    params - struct containing model parameters
%    options - struct containing algorithm parameters
%
% Outputs:
%    checks - struct
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: checkDynParameter

% Authors:       Maximilian Perschl
% Written:       17-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% search for checks in options 
switch errmsgid % TODO sort and categorize
    case 'change_intsteps'
        options = aux_changeOption_timeStep(params,options);
    otherwise
        CORAwarning('CORA:contDynamics','Unknown error ID. %s', errmsgid); return;
end


end


% Auxiliary functions -----------------------------------------------------

function options = aux_changeOption_timeStep(params,options)
    % round divisor up so time-step size only gets smaller
    % as more precise results are desired
    dividor = ceil(params.tFinal/options.timeStep);
    options.timeStep = params.tFinal/dividor;
    CORAwarning("CORA:global","options.timeStep has to divide params.tFinal in an integer number of steps." + ...
        " Therefore, options.timeStep was changed to %d",options.timeStep);
end

% ------------------------------ END OF CODE ------------------------------
