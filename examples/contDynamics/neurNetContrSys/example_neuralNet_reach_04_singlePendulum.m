function [completed,res,tTotal] = example_neuralNet_reach_04_singlePendulum
% example_neuralNet_reach_04_singlePendulum - example of reachability
%    analysis for a neural network controlled pendulum
%
% Syntax:
%    completed = example_neuralNet_reach_04_singlePendulum()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%    res - verification result
%    tTotal - total time
%
% Reference:
%   [1] Johnson, Taylor T., et al. "ARCH-COMP21 Category Report:
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for Continuous and Hybrid Systems Plants."
%       EPiC Series in Computing 80 (2021): 90-119.

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       08-November-2021
% Last update:   23-May-2022 (TL, ARCH'22 revisions)
%                30-March-2023 (TL, verify violated runs, ARCH'23 revisions)
% Last revision: 14-November-2022 (TL, clean up)

% ------------------------------ BEGIN CODE -------------------------------

disp("BENCHMARK: Single Pendulum")
rng(1)

% Parameters --------------------------------------------------------------

R0 = interval([1; 0], [1.2; 0.2]);

params.tFinal = 1;
params.R0 = polyZonotope(R0);

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.05;
options.alg = 'lin';
options.tensorOrder = 2;
options.taylorTerms = 4;
options.zonotopeOrder = 50;

% Parameters for NN evaluation --------------------------------------------

evParams = struct();
evParams.poly_method = "singh";

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
specUnsafe = specification(safeSet * 0.5 + 1, 'unsafeSet', spec.time);

% Verification ------------------------------------------------------------

t = tic;
[res, R, simRes] = verify(sys, spec, params, options, evParams, true);
tTotal = toc(t);
disp(['Result: ' res])

% Visualization -----------------------------------------------------------

disp("Plotting..")
figure; hold on; box on;

% plot specifications
plotOverTime(spec, 1, 'DisplayName', 'Safe set');
plotOverTime(specUnsafe, 1, 'DisplayName', 'Unsafe set');

% plot reachable set
useCORAcolors("CORA:contDynamics")
plotOverTime(R, 1, 'DisplayName', 'Reachable set');
updateColorIndex(); % don't plot initial set
% plotOverTime(R(1).R0, 1, 'DisplayName', 'Initial set');

% plot simulations
plotOverTime(simRes, 1, 'DisplayName', 'Simulations');

% labels and legend
xlabel('time');
ylabel('\theta');
legend();


% example completed -------------------------------------------------------

completed = true;

% handling for ARCH competition
if nargout < 2
    clear res;
end
if nargout < 3
    clear tTotal;
end

end

% ------------------------------ END OF CODE ------------------------------
