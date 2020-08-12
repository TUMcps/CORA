function [t,x,ind] = simulate(obj,params,options)
% simulate - simulates the system within a location
%
% Syntax:  
%    [t,x,ind] = simulate(obj,params,options)
%
% Inputs:
%    obj - nonlinParamSys object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .x0: initial point
%    options - ODE45 options (for hybrid systems)
%
% Outputs:
%    t - time vector
%    x - state vector
%    ind - returns the event which has been detected
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      26-May-2011 
% Last update:  08-May-2020 (MW, update interface)
% Last revision:---

%------------- BEGIN CODE --------------

if isempty(options.Events)
    [t,x] = ode45(getfcn(obj,params),...
        [params.tStart,params.tFinal],params.x0,options);
    ind = [];
else
    [t,x,te,xe,ind] = ode45(getfcn(obj,params),...
        [params.tStart,params.tFinal],params.x0,options);
end

%------------- END OF CODE --------------