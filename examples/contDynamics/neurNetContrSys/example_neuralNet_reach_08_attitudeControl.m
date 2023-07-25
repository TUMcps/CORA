function completed = example_neuralNet_reach_08_attitudeControl
% example_neuralNet_reach_08_attitudeControl - example of reachability
%    analysis for an neural network controlled system
%
% Syntax:
%    completed = example_neuralNet_reach_08_attitudeControl()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Reference:
%   [1] Johnson, Taylor T., et al. "ARCH-COMP22 Category Report:
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for Continuous and Hybrid Systems Plants."
%       EPiC Series in Computing TBD (2022): TBD.

% Author:       Tobias Ladner
% Written:      15-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

disp("BENCHMARK: Attitude Control")

% Parameters --------------------------------------------------------------

params.tFinal = 3;
params.R0 = polyZonotope(interval( ...
    [-0.45, -0.55, 0.65, -0.75, 0.85, -0.65], ...
    [-0.43, -0.53, 0.65, -0.74, 0.86, -0.64] ...
)');


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
f = @dynamics_attitudeControl;
sys = nonlinearSys(f);

% load neural network controller
% [4, 500, 2]
nn = neuralNetwork.readONNXNetwork('attitude_control_3_64_torch.onnx');

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.1);


% Specification -----------------------------------------------------------

unsafeSet = interval( ...
    [-0.2;-0.5;0;-0.7;0.7;-0.4], ...
    [0;-0.4;0.2;-0.6;0.8;-0.2] ...
);


% Simulation --------------------------------------------------------------

tic
simRes = simulateRandom(sys, params);
tSim = toc;
disp(['Time to compute random simulations: ', num2str(tSim)]);


% Check Violation --------------------------------------------------------

tic
isVio = false;
for i = 1:length(simRes)
    isVio = isVio || unsafeSet.contains(simRes(i).x{1}(end, :)');
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
    Rend = R(end).timePoint.set{end};
    Rend = interval(Rend);

    isVeri = ~isIntersecting(Rend, unsafeSet);
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

disp("Plotting..");
spec = specification(unsafeSet,'unsafeSet',interval(params.tFinal));

for w=1:3
    figure; hold on; box off;
    
    % plot specification 
    plotOverTime(spec,w, 'DisplayName', 'Unsafe set');

    % plot reachable set
    useCORAcolors('CORA:contDynamics')
    plotOverTime(R, w, 'DisplayName', 'Reachable set');

    % plot initial set
    plotOverTime(R(1).R0, w, 'DisplayName', 'Initial set');

    % plot simulation
    plotOverTime(simRes, w, 'DisplayName', 'Simulations');

    % labels and legend
    xlabel('time');
    ylabel(sprintf("\\omega_%d", w));
    legend(Location="best")
end

% example completed
completed = true;

%------------- END OF CODE --------------
