function [completed,res,tTotal] = example_neuralNet_reach_08_attitudeControl
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

disp("BENCHMARK: Attitude Control")

% Parameters --------------------------------------------------------------

params.tFinal = 3;
params.R0 = polyZonotope(interval( ...
    [-0.45, -0.55, 0.65, -0.75, 0.85, -0.65], ...
    [-0.44, -0.54, 0.66, -0.74, 0.86, -0.64] ...
)');

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.1;
options.alg = 'lin';
options.tensorOrder = 2;
options.taylorTerms = 4;
options.zonotopeOrder = 50;

% Parameters for NN evaluation --------------------------------------------

evParams = struct();
evParams.poly_method = "singh";

% System Dynamics ---------------------------------------------------------

% open-loop system
f = @dynamics_attitudeControl;
sys = nonlinearSys(f);

% load neural network controller
% [6, 64, 64, 64, 3]
nn = neuralNetwork.readONNXNetwork('attitude_control_3_64_torch.onnx');

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.1);


% Specification -----------------------------------------------------------

unsafeSet = interval( ...
    [-0.2;-0.5;0;-0.7;0.7;-0.4], ...
    [0;-0.4;0.2;-0.6;0.8;-0.2] ...
);
spec = specification(unsafeSet,'unsafeSet',interval(params.tFinal));

% Verification ------------------------------------------------------------

t = tic;
[res, R, simRes] = verify(sys, spec, params, options, evParams, true);
tTotal = toc(t);
disp(['Result: ' res])

% Visualization -----------------------------------------------------------

disp("Plotting..");

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
    legend(Location="north")
end


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
