function completed = example_neuralNet_reach_03_TORA
% example_neuralNet_reach_03_TORA - example of reachability analysis for a
%    neural network controlled cart
%
%
% Syntax:
%    completed = example_neuralNet_reach_03_TORA()
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

disp("BENCHMARK: Sherlock-Benchmark-9 (TORA)")

% Parameters --------------------------------------------------------------

R0 = interval([0.6;-0.7;-0.4;0.5],[0.7;-0.6;-0.3;0.6]);

params.tFinal = 20;
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
f = @(x,u) [x(2); -x(1) + 0.1*sin(x(3)); x(4); u(1) - 10];
sys = nonlinearSys(f);

% load neural network controller
% [4, 100, 100, 100, 1]
load('controllerTORA.mat');
nn = neuralNetworkOld(W, b, 'ReLU');

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 1);


% Specification -----------------------------------------------------------

safeSet = 2 * interval(-ones(4, 1), ones(4, 1));
spec = specification(safeSet, 'safeSet');


% Simulation --------------------------------------------------------------

tic
simRes = simulateRandom(sys, params);
tSim = toc;
disp(['Time to compute random simulations: ', num2str(tSim)]);


% Check Violation --------------------------------------------------------

tic
isVio = false;
for i = 1:length(simRes)
    x = simRes(i).x{1};
    for j=1:length(safeSet)
        isVio = isVio || ~all( ...
            (infimum(safeSet(j)) <= x(:, j)) & ...
            (x(:, j) <= supremum(safeSet(j))));
    end
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
    for i = 1:length(R)
        R_i = R(i);
        for j = 1:length(R_i.timeInterval)
            isVeri = isVeri & safeSet.contains(R_i.timeInterval.set{j});
        end
    end
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
xlim([-2.5, 2.5]);
ylim([-2.5, 2.5]);
legend()

figure; hold on; box on;
projDims = [3, 4];

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
xlabel('$x_3$ (angle)', 'interpreter', 'latex');
ylabel('$x_4\ (\dot{x_3})$', 'interpreter', 'latex')
xlim([-2.5, 2.5]);
ylim([-2.5, 2.5]);
legend()

% example completed
completed = true;

%------------- END OF CODE --------------