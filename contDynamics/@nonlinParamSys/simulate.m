function [t,x,ind,y] = simulate(obj,params,varargin)
% simulate - simulates the system within a location
%
% Syntax:  
%    [t,x] = simulate(obj,params)
%    [t,x,ind] = simulate(obj,params,options)
%    [t,x,ind,y] = simulate(obj,params,options)
%
% Inputs:
%    obj - nonlinParamSys object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .x0: initial point
%       .p: value for the parameter
%       .u: input signal
%    options - ODE45 options (for hybrid systems)
%
% Outputs:
%    t - time vector
%    x - state vector
%    ind - returns the event which has been detected
%    y - output vector
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

if nargout == 4
	warning("Output trajectories not supported for class nonlinParamSys!");
	y = [];
end

% parse input arguments
isOpt = false;
if nargin >= 3 && ~isempty(varargin{1})
   options = varargin{1}; 
   isOpt = true;
end

if ~isfield(params,'p')
    throw(CORAerror('CORA:specialError',...
        'Value for the parameter "params.p" is missing.'));
end

if ~isfield(params,'u')
	params.u = zeros(obj.nrOfInputs,1); 
end

if ~isfield(params,'tStart')
	params.tStart = 0;
end

tFinal = (params.tFinal-params.tStart)/size(params.u,2);

if isfield(params,'timeStep')
    tSpan = 0:params.timeStep:tFinal;
    if abs(tSpan(end)-tFinal) > 1e-10
       tSpan = [tSpan,tFinal]; 
    end
else
    tSpan = [0,tFinal];
end

% simulate the system
params_ = params;
t = [];
x = [];
ind = [];
x0 = params.x0;

for i = 1:size(params.u,2)

    params_.u = params.u(:,i);

    % simulate using MATLABs ode45 function
    try
        if isOpt
            [t_,x_,~,~,ind] = ode45(getfcn(obj,params_),tSpan,x0,options);
        else
            [t_,x_,~,~,ind] = ode45(getfcn(obj,params_),tSpan,x0);
        end
    catch
        if isOpt
            [t_,x_] = ode45(getfcn(obj,params_),tSpan,x0,options);
        else
            [t_,x_] = ode45(getfcn(obj,params_),tSpan,x0);
        end
    end

    % store the results
    x = [x;x_]; 
    if isempty(t)
        t = t_ + params.tStart;
    else 
        t = [t;t_ + t(end)]; 
    end
    x0 = x(end,:)';

    if ~isempty(ind)
       return; 
    end
end

%------------- END OF CODE --------------