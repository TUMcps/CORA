function completed = example_neuralNet_reach_02_ACC
% example_neuralNet_reach_02_ACC - example of reachability analysis for a
%    neural network adaptive cruise control
%
%
% Syntax:
%    completed = example_neuralNet_reach_02_ACC()
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
%   [2] Model Predictive Control Toolbox.
%       The MathWorks Inc., Natick, Massachusetts (2019).
%       https://www.mathworks.com/help/mpc/ug/adaptive-cruisecontrol-using-model-predictive-controller.html

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      08-November-2021
% Last update:  20-May-2022 (TL: ARCH'22 Revisions)
% Last revision:---

%------------- BEGIN CODE --------------

disp("BENCHMARK: Adaptive Cruise Controller (ACC)")

% Parameters --------------------------------------------------------------

R0 = interval([90; 32; 0; 10; 30; 0], [110; 32.2; 0; 11; 30.2; 0]);

params.tFinal = 5;
params.R0 = polyZonotope(R0);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.1;
options.taylorTerms = 4;
options.zonotopeOrder = 20;
options.alg = 'lin';
options.tensorOrder = 2;


% Parameters for NN evaluation --------------------------------------------
evParams = struct();
evParams.bound_approx = true;
evParams.polynomial_approx = "lin";


% System Dynamics ---------------------------------------------------------

% parameter
u_f = 0.0001;
a_lead = -2;
v_set = 30;
T_gap = 1.4;
D_default = 10;

% open-loop system
f = @(x, u) [x(2); x(3); -2 * x(3) + 2 * a_lead - u_f * x(2)^2; ...
    x(5); x(6); -2 * x(6) + 2 * u(1) - u_f * x(4)^2];
sys = nonlinearSys(f);

% affine map x_ = C*x + k mapping state x to input of neural network x_
C = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0; 1, 0, 0, -1, 0, 0; 0, 1, 0, 0, -1, 0];
k = [v_set; T_gap; 0; 0; 0];

% load neural network controller
% [5, 20, 20, 20, 20, 20, 1]
load('controllerACC.mat');

% distance computation
DM = [1 0 0 -1 0 0; 0 0 0 0 T_gap 0];
Db = [0;D_default];

actFun = [{'identity'}, repmat({'ReLU'}, [1, length(W)])];
W = [{C}, W];
b = [{k}, b];

nn = neuralNetworkOld(W, b, actFun);

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.1);


% Simulation --------------------------------------------------------------

tic;
simRes = simulateRandom(sys, params);
tSim = toc;
disp(['Time to compute random simulations: ', num2str(tSim)]);

% Check Violation --------------------------------------------------------

tic;
simResDistances = [];
isVio = false;
for i = 1:length(simRes)
    x = simRes(i).x{1};
    x = DM * x' + Db;
    
    simResDistances_i = simResult({x'}, simRes(i).t);
    if isempty(simResDistances)
        simResDistances = simResDistances_i;
    else
        simResDistances = add(simResDistances, simResDistances_i);
    end
        
    % relative distance D_rel
    distance = x(1, :);
    safe_distance = x(2, :);
    % safe distance D_safe
    isVio = isVio || ~all(distance >= safe_distance);
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

    tic;
    R = reach(sys, params, options, evParams);
    tComp = toc;
    disp(['Time to compute reachable set: ', num2str(tComp)]);

    % Verification --------------------------------------------------------

    tic;

    % transform 
    R_distances = DM * R + Db;

    isVeri = true;
    for i = 1:length(R_distances)
        for j = 1:length(R_distances(i).timeInterval.set)
            % read distances
            R_ij = R_distances(i).timeInterval.set{j};
            distance = interval(project(R_ij, 1));
            safe_distance = interval(project(R_ij, 2));

            % check safety
            isVeri = isVeri && (infimum(distance) > supremum(safe_distance));
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

% plot reachable sets

% relative distance D_rel
plotOverTime(R_distances, 1, 'DisplayName', 'Distance', 'FaceColor', CORAcolor("CORA:safe"));

% safe distance D_safe
plotOverTime(R_distances, 2, 'DisplayName', 'Safe distance', 'FaceColor', CORAcolor("CORA:unsafe"));

% plot simulations
plotOverTime(simResDistances, 1, 'k', 'DisplayName', 'Simulations');

% labels and legend
xlabel('time');
ylabel('distance');
legend();

% example completed
completed = true;

%------------- END OF CODE --------------
