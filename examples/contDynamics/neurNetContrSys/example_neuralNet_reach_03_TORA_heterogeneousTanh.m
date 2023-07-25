function completed = example_neuralNet_reach_03_TORA_heterogeneousTanh
% example_neuralNet_reach_03_TORA_heterogeneousTanh - example of
%    reachability analysis for a neural network controlled cart
%
% Syntax:
%    completed = example_neuralNet_reach_03_TORA_heterogeneousTanh()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Reference:
%   [1] Johnson, Taylor T., et al. "ARCH-COMP21 Category Report:
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for Continuous and Hybrid Systems Plants."
%       EPiC Series in Computing 80 (2021): 90-119.
%   [2] Dutta, Souradeep, et al. "Reachability analysis for
%       neural feedback systems using regressive polynomial rule inference"
%       Proceedings of the 22nd ACM International Conference on
%       Hybrid Systems: Computation and Control. 2019.
%   [3] Jankovic, Mrdjan, et al. "TORA example:
%       cascade-and passivity-based control designs."
%       IEEE Transactions on Control Systems Technology 4.3 (1996): 292-297

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      08-November-2021
% Last update:  20-May-2022 (TL: ARCH'22 Revisions)
% Last revision:---

%------------- BEGIN CODE --------------

disp("BENCHMARK: TORA Heterogeneous (tanh)")

% Parameters --------------------------------------------------------------

R0 = interval([-0.77;-0.45;0.51;-0.3],[-0.75;-0.43;0.54;-0.28]);

params.tFinal = 5;
params.R0 = polyZonotope(R0);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.05;
options.taylorTerms = 4;
options.zonotopeOrder = 200;
options.alg = 'lin';
options.tensorOrder = 3;
options.errorOrder = 10;
options.intermediateOrder = 50;


% Parameters for NN evaluation --------------------------------------------

evParams = struct();
evParams.bound_approx = true;
evParams.polynomial_approx = "lin";


% System Dynamics ---------------------------------------------------------

% open-loop system
f = @(x,u) [x(2); -x(1) + 0.1*sin(x(3)); x(4); 11*u(1)];
sys = nonlinearSys(f);

% load neural network controller
% [4, 100, 100, 100, 1]
load('nn_tora_relu_tanh.mat');
nn = neuralNetworkOld(W, b, {'ReLU', 'ReLU', 'ReLU', 'tanh'});

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.5);


% Specification -----------------------------------------------------------

goalSet = interval([-0.1;-0.9;-Inf;-Inf], [0.2;-0.6;Inf;Inf]);
spec = specification(goalSet, 'safeSet');


% Simulation --------------------------------------------------------------

tic
simRes = simulateRandom(sys, params);
tSim = toc;
disp(['Time to compute random simulations: ', num2str(tSim)]);

% Check Violation --------------------------------------------------------

tic
isVio = false;
for i = 1:length(simRes)
    isVio = isVio || ~goalSet.contains(simRes(i).x{1}(end, :)');
end
tVio = toc;
disp(['Time to check violation in simulations: ', num2str(tVio)]);

if isVio
    disp("Result: VIOLATED")
    R = params.R0;
    tComp = 0;
    tVeri = 0;
else
    % Reachability Analysis -----------------------------------------------

    tic
    R = reach(sys, params, options, evParams);
    tComp = toc;
    disp(['Time to compute reachable set: ', num2str(tComp)]);

    % Verification --------------------------------------------------------

    tic
    isVeri = true;
    isVeri = isVeri & goalSet.contains(R(end).timePoint.set{end});
    tVeri = toc;
    disp(['Time to check verification: ', num2str(tVeri)]);

    if isVeri
        disp('Result: VERIFIED');
    else
        disp('Result: UNKNOWN');
    end
end
disp(['Total Time: ', num2str(tSim+tVio+tComp+tVeri)]);


% Visualization -----------------------------------------------------------

disp("Plotting..")
figure; hold on; box on;
projDims = [1, 2];

% plot specifications
plot(spec, projDims, 'DisplayName', 'Safe set');

% plot reachable set
useCORAcolors("CORA:contDynamics")
 plot(R, projDims, 'DisplayName', 'Reachable set');

% plot initial set
plot(R(1).R0, projDims, 'DisplayName', 'Initial set');

% plot simulations
plot(simRes, projDims, 'DisplayName', 'Simulations');

% labels, limits, and legend
xlabel('$x_1$ (distance)', 'interpreter', 'latex');
ylabel('$x_2\ (\dot{x_1})$', 'interpreter', 'latex')
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
legend()

% example completed
completed = true;

%------------- END OF CODE --------------