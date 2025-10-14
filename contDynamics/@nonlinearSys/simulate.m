function [t,x,ind,y] = simulate(nlnsys,params,varargin)
% simulate - simulates the system within a location
%
% Syntax:
%    [t,x] = simulate(nlnsys,params)
%    [t,x,ind] = simulate(nlnsys,params,options)
%    [t,x,ind,y] = simulate(nlnsys,params,options)
%
% Inputs:
%    nlnsys - nonlinearSys object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .timeStep: time step size
%       .x0: initial point
%       .u: input
%    options - ODE45 options (for hybrid systems)
%
% Outputs:
%    t - time vector
%    x - state vector
%    ind - returns the event which has been detected
%    y - output vector
%
% Example: 
%    nlnsys = nonlinearSys(@vanderPolEq);
%
%    params.x0 = [1.4;2.3];
%    params.tFinal = 6;
%
%    [t,x] = simulate(nlnsys,params);
%
%    plot(x(1,:),x(2,:),'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       03-May-2007 
% Last update:   12-March-2008
%                08-May-2020 (MW, update interface)
%                28-August-2025 (LL, transpose t, x, and y)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargout == 4
    CORAwarning('CORA:contDynamics',"Output trajectories not supported for class nonlinearSys!");
	y = [];
end

% parse input arguments
isOpt = false;
if nargin >= 3 && ~isempty(varargin{1})
	options = varargin{1}; 
    isOpt = true;
end

% default initialization of input by an all-zero vector 
if ~isfield(params,'u')
	params.u = zeros(nlnsys.nrOfInputs,1); 
end

% default initialization of time by 0
if ~isfield(params,'tStart')
	params.tStart = 0;
end

% tFinal is evenly split in case the input changes, so that each input is
% active for the same time
tFinal = (params.tFinal-params.tStart)/size(params.u,2);

% determine time span for simulation
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
    % simulate until a zero-crossing is detected
    try
        % options are specified
        if isOpt
            [t_,x_,~,~,ind] = ode45(getfcn(nlnsys,params_),tSpan,x0,options);
        % options are not specified
        else
            [t_,x_,~,~,ind] = ode45(getfcn(nlnsys,params_),tSpan,x0);
        end
    % no zero-crossing specified
    catch
        % options are specified
        if isOpt
            [t_,x_] = ode45(getfcn(nlnsys,params_),tSpan,x0,options);
        % options are not specified
        else
            [t_,x_] = ode45(getfcn(nlnsys,params_),tSpan,x0);
        end
    end

    % store the results
    x = [x x_']; 
    if isempty(t)
        t = t_' + params.tStart;
    else 
        t = [t t_'+t(end)]; 
    end
    x0 = x(:,end);

    if ~isempty(ind)
        return; 
    end
end

% ------------------------------ END OF CODE ------------------------------
