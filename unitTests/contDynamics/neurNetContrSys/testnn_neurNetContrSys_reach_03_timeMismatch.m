function res = testnn_neurNetContrSys_reach_03_timeMismatch
% testnn_neurNetContrSys_reach_03_timeMismatch - unit test function of 
%    neurNetContrSys where reachability analysis time step and
%    nn sampling time don't match
%
% toy example for neurNetContrSys - car example
%
% Syntax:
%    res = testnn_neurNetContrSys_reach_03_timeMismatch
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Tobias Ladner
% Written:       20-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 1;
% [location; velocity; acceleration]
params.R0 = polyZonotope(interval([100;50;0], [110;60;0]));
params.U = interval(0, 0);

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.04;
options.taylorTerms = 4;
options.zonotopeOrder = 50;
options.alg = 'lin';
options.tensorOrder = 2;

% Parameters for NN evaluation --------------------------------------------
evParams = struct();
evParams.bound_approx = true;
evParams.poly_method = "singh";

% System Dynamics ---------------------------------------------------------

% open-loop system
% Simple Car dynamics
dynamics = @(x, u) [; ...
    x(2); ... % \dot{l} = v
    x(3); ... % \dot{v} = a
    u(1) ... % \dot{a} = u | input of driver
    ];
sys = nonlinearSys(dynamics);

% load neural network controller
neurons = [3, 10, 1];

% random weights
W = cell(length(neurons)-1, 1);
b = cell(length(neurons)-1, 1);
scale = 4;
for i=1:length(neurons)-1
    W{i} = rand(neurons(i+1), neurons(i)) * scale - scale/2;
    b{i} = rand(neurons(i+1), 1) * scale - scale/2;
end

% init layers
layers = cell(length(W)*2, 1);
for i=1:length(W)
    layers{2*i-1} = nnLinearLayer(W{i}, b{i});
    layers{2*i} = nnSigmoidLayer();
end
nn = neuralNetwork(layers);

% construct neural network controlled system
dt = 0.1; % nn sampling time
sys = neurNetContrSys(sys, nn, dt);

% Run ---------------------------------------------------------

% compute reachable set
R = reach(sys, params, options, evParams);

% check if subsequent time intervals have matching time bounds
T = [R.timeInterval];
T = [T.time];

resvec = [];
resvec(end+1) = T{1}.inf == 0;
for i=1:length(T)-1
    resvec(end+1) = T{i}.sup == T{i+1}.inf;
end
resvec(end+1) = T{end}.sup == params.tFinal;

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
