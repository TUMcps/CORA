function [t, x, ind, y] = simulate(obj, params, varargin)
% simulate - simulates a neural network controlled system
%
% Syntax:
%    [t,x] = simulate(obj,params)
%    [t,x,ind] = simulate(obj,params,options)
%    [t,x,ind,y] = simulate(obj,params,options)
%
% Inputs:
%    obj - neurNetContrSys object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time (default 0)
%       .tFinal: final time
%       .x0: initial point
%    options - ODE45 options (for hybrid systems)
%
% Outputs:
%    t - time vector
%    x - state vector
%    ind - returns the event which has been detected
%    y - output vector
%
% Example:
%    % dynamic system
%    f = @(x,u) [x(2) + u(2); (1-x(1)^2)*x(2) - x(1) + u(1)];
%    sysOL = nonlinearSys(f);
% 
%    % neural network controller
%    layers = cell(4, 1);
%    W1 = rand(10,2); b1 = rand(10,1);
%    layers{1} = nnLinearLayer(W1, b1);
%    layers{2} = nnSigmoidLayer();
%    W2 = rand(2,10); b2 = rand(2,1);
%    layers{3} = nnLinearLayer(W2, b2);
%    layers{4} = nnSigmoidLayer();
%    nn = neuralNetwork(layers);
% 
%    % neural network controlled system
%    dt = 0.01;
%    sys = neurNetContrSys(sysOL,nn,dt);
% 
%    params.x0 = [1;2];
%    params.tFinal = 1;
% 
%    [t,x] = simulate(sys,params);
% 
%    figure
%    plot(x(:,1),x(:,2),'k', 'DisplayName', 'Simulation');
%    legend()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neurNetContrSys, contDynamics/simulateRandom

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       13-December-2021
% Last update:   28-November-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargout == 3
    ind = [];
end
if nargout == 4
    warning("Output trajectories not supported for class nonlinearSysDT!");
    y = [];
end

isOpt = false;
if nargin >= 3 && ~isempty(varargin{1})
    options = varargin{1};
    isOpt = true;
end

if ~isfield(params, 'tStart')
    params.tStart = 0;
end

if ~isfield(params, 'u')
    params.u = obj.nn.evaluate(params.x0);
end

time = (params.tStart:obj.dt:params.tFinal)';
if time(end) ~= params.tFinal
    time = [time; params.tFinal];
end

% consider changing inputs
if size(params.u, 2) ~= 1
    tu = linspace(params.tStart, params.tFinal, size(params.u, 2));
else
    tu = time;
    params.u = params.u * ones(1, length(time)-1);
end

% simulate the system
params_ = params;
t = [];
x = [];
ind = [];
x0 = params.x0;
cnt = 1;

for i = 1:length(time) - 1

    % compute control input for the current state
    u = evaluate(obj.nn, x0);
    x0 = [x0; u];

    % get number of disturbance changes duing sampling time
    index = find(tu(cnt:end) < time(i+1));
    w = params.u(cnt:cnt+index(end)-1);

    % loop over the number of disturbance changes
    for j = 1:size(w, 2)

        params_.u = params.u(:, j);
        params_.x0 = x0;
        params_.w = zeros(length(x0), 1); % for linear systems

        if isfield(params, 'timeStep')
            tSpan = tu(cnt+j-1):params.timeStep:tu(cnt+j);
            if abs(tSpan(end)-tu(cnt+j)) > 1e-10
                tSpan = [tSpan, tu(cnt+j)];
            end
        else
            tSpan = [tu(cnt+j-1), tu(cnt+j)];
        end

        % simulate using MATLABs ode45 function
        try
            if isOpt
                [t_, x_, ~, ~, ind] = ode45(getfcn(obj.sys, params_), tSpan, ...
                    x0, options);
            else
                [t_, x_, ~, ~, ind] = ode45(getfcn(obj.sys, params_), tSpan, x0);
            end
        catch
            if isOpt
                [t_, x_] = ode45(getfcn(obj.sys, params_), tSpan, x0, options);
            else
                [t_, x_] = ode45(getfcn(obj.sys, params_), tSpan, x0);
            end
        end

        % store the results
        x = [x; x_(:, 1:obj.dim)];
        t = [t; t_];
        x0 = x(end, 1:obj.dim)';
        cnt = cnt + index(end);

        if ~isempty(ind)
            return;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
