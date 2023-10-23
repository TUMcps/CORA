function res = test_neurNetContrSys_reach_01_toyExample
% test_neurNetContrSys_reach_01_toyExample - unit test function of 
%    neurNetContrSys
%
% toy example for neurNetContrSys - car example
%
% Syntax:
%    res = test_neurNetContrSys_reach_01_toyExample
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Tobias Ladner
% Written:       02-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 1;
% [location; velocity; acceleration]
params.R0 = polyZonotope(interval([100;50;0], [110;60;0]));
params.U = interval(0, 0);

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.1;
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

% compute simulations
simRes = simulateRandom(sys, params);

% compute reachable set
R = reach(sys, params, options, evParams);

% plot
% figure; hold on;
% rs = plot(R, [1, 2], 'FaceColor', [.8, .8, .8], 'EdgeColor', 'none');
% is = plot(params.R0, [1, 2], 'FaceColor', 'w', 'EdgeColor', 'k');
% ss = plot(simRes, [1, 2], 'k');

res = true;

% ------------------------------ END OF CODE ------------------------------
