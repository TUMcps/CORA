function [t,x,ind,y] = simulate(obj,params,varargin)
% simulate - simulates the linear interval system within a location
%
% Syntax:
%    [t,x] = simulate(obj,params)
%    [t,x,ind] = simulate(obj,params,options)
%    [t,x,ind,y] = simulate(obj,params,options)
%
% Inputs:
%    obj - linParamSys object
%    params - struct containing the parameters for the simulation
%       .tStart initial time t0
%       .tFinal final time tf
%       .x0 initial point x0
%       .u piecewise constant input signal u(t) specified as a matrix
%           for which the number of rows is identical to the number of
%           system input
%    options - ODE45 options (for hybrid systems)
%
% Outputs:
%    t - time vector
%    x - state vector
%    ind - returns the event which has been detected
%    y - output vector
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       16-May-2007 
% Last update:   07-January-2009
%                08-May-2020 (MW, update interface)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargout == 4
    CORAwarning('CORA:contDynamics',"Output trajectories not supported for class linParamSys!");
    y = [];
end

% parse input arguments
isOpt = 0;
if nargin >= 3 && ~isempty(varargin{1})
	options = varargin{1}; 
	isOpt = true;
end

% default initialization of input by an all-zero vector 
if ~isfield(params,'u')
	params.u = zeros(obj.nrOfInputs,1); 
end

% default initialization of time by 0
if ~isfield(params,'tStart')
	params.tStart = 0;
end

% tFinal is evenly split in case the input changes, so that each input is
% active for the same time
tFinal = params.tFinal/size(params.u,2);

% determine time span for simulation
if isfield(params,'timeStep')
    tSpan = params.tStart:params.timeStep:tFinal;
    if abs(tSpan(end)-tFinal) > 1e-10
    	tSpan = [tSpan,tFinal]; 
    end
else
    tSpan = [params.tStart,params.tFinal];
end

% sample one system matrix
tmp = randPoint(obj.A);
obj.sampleMatrix.A = tmp;

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
            [t_,x_,~,~,ind] = ode45(getfcn(obj,params_),tSpan,x0,options);
        % options are not specified
        else
            [t_,x_,~,~,ind] = ode45(getfcn(obj,params_),tSpan,x0);
        end
    % no zero-crossing specified
    catch
        % options are specified
        if isOpt
            [t_,x_] = ode45(getfcn(obj,params_),tSpan,x0,options);
        % options are not specified
        else
            [t_,x_] = ode45(getfcn(obj,params_),tSpan,x0);
        end
    end

    % store the results
    x = [x;x_]; 
    if isempty(t)
        t = [t;t_];
    else
        t = [t;t_ + t(end)]; 
    end
    x0 = x(end,:)';

    if ~isempty(ind)
    	return; 
    end
end

% ------------------------------ END OF CODE ------------------------------
