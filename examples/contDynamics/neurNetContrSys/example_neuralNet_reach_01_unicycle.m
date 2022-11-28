function completed = example_neuralNet_reach_01_unicycle
% example_neuralNet_reach_01_unicycle - example of reachability analysis
%    for a neural network controlled system
%
% Syntax:
%    completed = example_neuralNet_reach_01_unicycle()
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
% Written:      17-September-2021
% Last update:  20-May-2022 (TL: ARCH'22 Revisions)
% Last revision:---

%------------- BEGIN CODE --------------

disp("BENCHMARK: Sherlock-Benchmark 10 (Unicycle Car Model)");

% Parameters --------------------------------------------------------------

params.tFinal = 10;
params.R0 = polyZonotope(interval([9.5; -4.5; 2.1; 1.5], [9.55; -4.45; 2.11; 1.51]));


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.1;
options.taylorTerms = 4;
options.zonotopeOrder = 50;
options.alg = 'lin';
options.tensorOrder = 2;


% Parameters for NN evaluation --------------------------------------------

evParams = struct();
evParams.bound_approx = true;
evParams.polynomial_approx = "lin";


% System Dynamics ---------------------------------------------------------

% open-loop system
f = @(x, u) [x(4) * cos(x(3)); x(4) * sin(x(3)); u(2) - 20; u(1) - 20 + u(3)];
sys = nonlinearSys(f);

% load neural network controller
% [4, 500, 2]
load('controllerUnicycle.mat');
nn = neuralNetwork({ ...
    nnLinearLayer(W{1}, b{1}), ...
    nnReLULayer(), ...
    nnLinearLayer(W{2}, b{2}), ...
    nnReLULayer(), ...
});

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.2);


% Specification -----------------------------------------------------------

goalSet = interval(-[0.6; 0.2; 0.06; 0.3], [0.6; 0.2; 0.06; 0.3]);


% Simulation --------------------------------------------------------------

tic
simRes = simulateRandom(sys, params);
tSim = toc;
disp(['Time to compute random simulations: ', num2str(tSim)]);


% Check Violation --------------------------------------------------------

tic
isVio = false;
for i = 1:length(simRes.x)
    isVio = isVio || ~goalSet.contains(simRes.x{i}(end, :)');
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

    tic;
    isVeri = goalSet.contains(R(end).timePoint.set{end});
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

% plot goal set
gs = plot(goalSet, [1, 2], 'FaceColor', [0, .8, 0]);

% plot reachable set
rs = plot(R, [1, 2], 'FaceColor', [.8, .8, .8]);

% plot initial set
is = plot(params.R0, [1, 2], 'k', 'FaceColor', 'w');

% plot simulations
ss = plot(simRes,[1, 2], 'k');

% labels and legend
xlabel('x'); ylabel('y');
legend([gs, rs, is, ss], "Goal Set", "Reachable Set", ...
    "Initial Set", "Simulations", Location="best")


% example completed
completed = true;

%------------- END OF CODE --------------

