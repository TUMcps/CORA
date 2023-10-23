function [completed,res,tTotal] = benchmark_neuralNet_reach_09_QUAD
% benchmark_neuralNet_reach_09_QUAD - example of reachability analysis
%    for an neural network controlled system
%
% Syntax:
%    completed = benchmark_neuralNet_reach_09_QUAD()
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean
%    res - verification result
%    tTotal - total time
%
% Reference:
%   [1] Lopez, Diego Manzanas, et al. "ARCH-COMP22 category report: 
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for continuous and hybrid systems plants." Proceedings of 
%       9th International Workshop on Applied Verification of 
%       Continuous and Hybrid Systems (ARCH22). 2022.

% Authors:       Tobias Ladner
% Written:       15-June-2022
% Last update:   30-March-2022 (TL, ARCH'23 revisions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp("BENCHMARK: Quadrotor (QUAD)")

% Parameter ---------------------------------------------------------------

params.tFinal = 5;
w = 0.4;
params.R0 = polyZonotope(interval( ...
    [-w; -w; -w; -w; -w; -w; 0; 0; 0; 0; 0; 0], ...
    [w; w; w; w; w; w; 0; 0; 0; 0; 0; 0] ...
));

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.01;
options.alg = 'lin';
options.tensorOrder = 2;
options.taylorTerms = 1;
options.zonotopeOrder = 80;

% Parameters for NN evaluation --------------------------------------------

evParams = struct();
evParams.poly_method = 'regression';

% System Dynamics ---------------------------------------------------------

% open-loop system
f = @dynamics_quad;
sys = nonlinearSys(f);

% load neural network controller
% [12, 64, 64, 64, 3]
nn = neuralNetwork.readONNXNetwork('quad_controller_3_64_torch.onnx');
nn.evaluate(params.R0, evParams);
nn.refine(2, "layer", "both", params.R0.c, true);

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.1);

% Specification -----------------------------------------------------------

goalSet = interval( ...
    [-Inf;-Inf; 0.94;-Inf;-Inf;-Inf;-Inf;-Inf;-Inf;-Inf;-Inf;-Inf], ...
    [ Inf; Inf; 1.06; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf] ...
);
spec = specification(goalSet, 'safeSet', interval(params.tFinal));

% Verification ------------------------------------------------------------

t = tic;
[res, R, simRes] = verify(sys, spec, params, options, evParams, true);
tTotal = toc(t);
disp(['Result: ' res])

% Visualization -----------------------------------------------------------
disp("Plotting..")

figure; hold on; box on;

% plot specification (over entire time horizon)
spec = specification(goalSet, 'safeSet', interval(0, params.tFinal));
plotOverTime(spec, 3, 'DisplayName', 'Goal set');

% plot reachable set
useCORAcolors('CORA:contDynamics')
plotOverTime(R, 3, 'DisplayName', 'Reachable set');

% plot initial set
plotOverTime(R(1).R0, 3, 'DisplayName', 'Initial set');

% plot simulation
plotOverTime(simRes, 3, 'DisplayName', 'Simulations');

% labels and legend
xlabel('time'); ylabel('altitude');
ylim([-0.5, 2])
legend('Location', "southeast")


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
