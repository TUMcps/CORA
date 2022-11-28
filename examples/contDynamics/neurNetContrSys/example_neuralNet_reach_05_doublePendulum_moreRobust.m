function completed = example_neuralNet_reach_05_doublePendulum_moreRobust
% example_neuralNet_reach_05_doublePendulum_moreRobust - example of
%    reachability analysis for a neural network controlled double pendulum
%                                  
%
% Syntax:  
%    completed = example_neuralNet_reach_05_doublePendulum_moreRobust()
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

disp("BENCHMARK: Double Pendulum (more robust)")

% Parameters --------------------------------------------------------------

R0 = interval([1;1;1;1],[1.3;1.3;1.3;1.3]);

params.tFinal = 0.4;
params.R0 = polyZonotope(R0);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.02;
options.taylorTerms = 4;
options.zonotopeOrder = 200;
options.alg = 'lin';
options.tensorOrder = 2;
options.lagrangeRem.simplify = 'optimize';


% Parameters for NN evaluation --------------------------------------------

evParams = struct();
evParams.bound_approx = true;
evParams.polynomial_approx = "lin";


% System Dynamics ---------------------------------------------------------

% open-loop system
sys = nonlinearSys(@doublePendulum);

% load neural network controller
% [4, 25, 25, 2]
nn = neuralNetwork.readONNXNetwork('controller_double_pendulum_more_robust.onnx');

% construct neural network controlled system
sys = neurNetContrSys(sys,nn,0.02);


% Specification -----------------------------------------------------------

safeSet = interval([-0.5;-0.5;-0.5;-0.5],[1.5;1.5;1.5;1.5]);
spec = specification(safeSet,'safeSet',interval(0.5,1));


% Simulation --------------------------------------------------------------

tic
simRes = simulateRandom(sys, params);
tSim = toc;
disp(['Time to compute random simulations: ', num2str(tSim)]);


% Check Violation --------------------------------------------------------

tic
isVio = false;
for i = 1:length(simRes.x)
    x = simRes.x{i};
    for j =1:length(safeSet)
        isVio = isVio || ~all( ...
            (infimum(safeSet(j)) <= x(:, j)) & ...
            (x(:, j) <= supremum(safeSet(j))));
    end
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
    disp(['Time to check Verification: ', num2str(tVeri)]);

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

% plot specification
ss = plot(safeSet, [1,2], 'FaceColor', [0,.8,0]);

% plot reachable set
if ~isVio
    plot(R, [1, 2],'FaceColor',[0.7 0.7 0.7]);
end

% plot initial set
is = plot(R0,[1,2],'FaceColor','w','EdgeColor','k');

% plot simulations
sims = plot(simRes,[1,2],'k');

% labels and legend
xlabel('\theta_1'); ylabel('\theta_2');
legend([ss, is, sims], "Safe Set", "Initial Set", "Simulations", Location="northwest")


figure; hold on; box on;

% plot specifications
ss = plot(safeSet, [3,4], 'FaceColor',[0,.8,0]);

% plot reachable set
if ~isVio
    plot(R, [3, 4],'FaceColor',[0.7 0.7 0.7]);
end

% plot initial set
is = plot(R0,[3,4],'FaceColor','w','EdgeColor','k');

% plot simulations
sims = plot(simRes,[3,4],'k');

% labels and legend
xlabel('$\dot \theta_1$','interpreter','latex');
ylabel('$\dot \theta_2$','interpreter','latex');
legend([ss, is, sims], "Safe Set", "Initial Set", "Simulations", Location="northwest")

% example completed
completed = true;

%------------- END OF CODE --------------