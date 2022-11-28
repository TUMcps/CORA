function completed = example_neuralNet_reach_04_singlePendulum
% example_neuralNet_reach_04_singlePendulum - example of reachability
%    analysis for a neural network controlled pendulum
%
%
% Syntax:
%    completed = example_neuralNet_reach_04_singlePendulum()
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

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      08-November-2021
% Last update:  23-May-2022 (TL: ARCH'22 Revisions)
% Last revision:14-November-2022 (TL: clean up)

%------------- BEGIN CODE --------------

disp("BENCHMARK: Single Pendulum")

% Parameters --------------------------------------------------------------

R0 = interval([1; 0], [1.2; 0.2]);

params.tFinal = 1;
params.R0 = polyZonotope(R0);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.05;
options.taylorTerms = 4;
options.zonotopeOrder = 50;
options.alg = 'lin';
options.tensorOrder = 2;


% Parameters for NN evaluation --------------------------------------------

evParams = struct();
evParams.bound_approx = true;
evParams.polynomial_approx = "lin";


% System Dynamics ---------------------------------------------------------

% parameters
m = 0.5;
L = 0.5;
c = 0;
g = 1;

% open-loop system (u = T)
f = @(x, u) [x(2); g / L * sin(x(1)) + 1 / (m * L^2) * (u(1) - c * x(2))];
sys = nonlinearSys(f);

% load neural network controller
% [2, 25, 25, 1]
nn = neuralNetwork.readONNXNetwork('controller_single_pendulum.onnx');

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.05);


% Specification -----------------------------------------------------------

safeSet = interval([0; -Inf], [1; Inf]);
spec = specification(safeSet, 'safeSet', interval(0.5, 1));


% Simulation --------------------------------------------------------------

tic;
simRes = simulateRandom(sys, params);
tSim = toc;
disp(['Time to compute random simulations: ', num2str(tSim)]);


% Check Violation --------------------------------------------------------

tic
isVio = false;
for i = 1:length(simRes.x)
    x = simRes.x{i}(simRes.t{i} >= 0.5);
    isVio = isVio || ~all((0 <= x(:, 1)) & (x(:, 1) <= 1));
end
tVio = toc;
disp(['Time to check violation in simulations: ', num2str(tVio)]);

if isVio
    disp("Result: VIOLATED")
    R = [];
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

% plot specifications
ss = plot(cartProd(spec.time, interval(spec.set)), [1, 2], 'FaceColor', [0, .8, 0]);
us = plot(cartProd(spec.time, 1 + 0.4*interval(spec.set)), [1, 2], 'FaceColor', [.8, 0, 0]);
alpha(us,.5)

% plot reachable set
if ~isVio
    plotOverTime(R, 1, 'FaceColor', [0.7, 0.7, 0.7]);
end

% plot simulations
plotOverTime(simRes, 1, 'k');

% labels and legend
xlabel('time');
ylabel('\theta');
legend([us, ss], "Unsafe Set", "Safe Set", Location="best");

% example completed
completed = true;

%------------- END OF CODE --------------