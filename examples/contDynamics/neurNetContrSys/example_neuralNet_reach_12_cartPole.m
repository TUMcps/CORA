function [completed,res,tTotal] = example_neuralNet_reach_12_cartPole
% example_neuralNet_reach_12_cartPole - example of reachability analysis
%    for a neural network controlled system
%
% Syntax:
%    completed = example_neuralNet_reach_12_cartPole()
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
%   [1] "ARCH-COMP24 Category Report:
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for Continuous and Hybrid Systems Plants."

% Authors:       Tobias Ladner
% Written:       24-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp("BENCHMARK: Cartpole");

% Parameters --------------------------------------------------------------

params.tFinal = 10;
params.R0 = polyZonotope( ...
    enlarge(interval([-0.1; -0.05; -0.1; -0.05], [0.1; 0.05; 0.1; 0.05]),0.4) ...
);

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.02;
options.alg = 'lin';
options.tensorOrder = 2;
options.taylorTerms = 1;
% options.errorOrder = 10;
% options.intermediateOrder = 10;
options.zonotopeOrder = 50;

% Options for NN evaluation -----------------------------------------------

options.nn = struct();
options.nn.poly_method = "regression";

% System Dynamics ---------------------------------------------------------

% open-loop system
f = @(x, u) [x(2); 2*u(1); x(4); (0.08*0.41*(9.8*sin(x(3))-2*u(1)*cos(x(3)))-0.0021*x(4))/0.0105];
sys = nonlinearSys(f);

% load neural network controller
% [4, 64, 64, 1]
nn = neuralNetwork.readONNXNetwork('model-cartPole.onnx');
% nn.evaluate(params.R0);
% nn.refine();

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.02);

% Specification -----------------------------------------------------------

goalSet = interval(-[0.001; inf; 0.001; 0.001], [0.001; inf; 0.001; 0.001]);
spec = specification(goalSet, 'safeSet', interval(8,params.tFinal));

% Verification ------------------------------------------------------------

t = tic;
[res, R, traj] = verify(sys, spec, params, options, true, 10);
tTotal = toc(t);
disp(['Result: ' res])

% Visualization -----------------------------------------------------------

disp("Plotting..")
figure; hold on; box on;

% plot reachable set
useCORAcolors("CORA:contDynamics")
plot(R, [1, 3], 'DisplayName', 'Reachable set');

% plot initial set
plot(R(1).R0, [1, 3], 'DisplayName', 'Initial set');

% plot goal set
plot(specification(goalSet, 'safeSet'), [1, 3], 'DisplayName', 'Goal set');
updateColorIndex(2)

% plot simulations
plot(traj,[1, 3], 'DisplayName', 'Simulations');

% labels and legend
xlabel('x_{(1)} (position)'); ylabel('x_{(3)} (angle)');
legend()

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
