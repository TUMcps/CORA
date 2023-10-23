function [completed,res,tTotal] = example_neuralNet_reach_10_spacecraftDocking
% example_neuralNet_reach_10_spacecraftDocking - example of reachability analysis
%    for an neural network controlled system
%
% Syntax:
%    completed = example_neuralNet_reach_10_spacecraftDocking()
%
% Inputs:
%    no
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
% Written:       20-June-2022
% Last update:   30-March-2022 (TL, ARCH'23 revisions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp("BENCHMARK: Spacecraft Docking")

% Parameter ---------------------------------------------------------------

params.tFinal = 40;
params.R0 = polyZonotope(interval( ...
    [70; 70; -0.28; -0.28], ...
    [106; 106; 0.28; 0.28] ...
));

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.1;
options.alg = 'lin';
options.tensorOrder = 3;
options.taylorTerms = 4;
options.zonotopeOrder = 200;
options.errorOrder = 10;
options.intermediateOrder = 50;

% Parameters for NN evaluation --------------------------------------------

evParams = struct();
evParams.poly_method = "singh";

% System Dynamics ---------------------------------------------------------

% open-loop system
f = @dynamics_spacecraftDocking;
sys = nonlinearSys(f);

% load neural network controller
% [4, 256, 256, 2]
nn = neuralNetwork.readONNXNetwork('controller_spacecraftDocking.onnx');
% nn.evaluate(params.R0, evParams);
% nn.refine(2, "layer", "both", params.R0.c, true);


% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.1);

% Specification -----------------------------------------------------------

v0 = 0.2;
v1 = 2*0.001027;
isSafe = @(x) (sqrt(x(3)^2 + x(4)^2)) <= (v0 + v1*sqrt(x(1)^2 + x(2)^2));

% Simulation --------------------------------------------------------------

tic
simRes = simulateRandom(sys, params);
tSim = toc;
disp(['Time to compute random simulations: ', num2str(tSim)]);

% Check Violation --------------------------------------------------------

tic
isVio = false;
for i = 1:length(simRes)
    x_i = simRes(i).x;
    for j = 1:length(x_i)
        x_ij = x_i{j};
        for k = size(x_ij,1)
            isVio = isVio || ~isSafe(x_ij(k, :)');
            
            if isVio
                break;
            end
        end
    end
end
tVio = toc;
disp(['Time to check violation in simulations: ', num2str(tVio)]);


if isVio 
    res = 'VIOLATED';
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
    Rend = R(end).timePoint.set{end};
    % Rend = interval(Rend);
    %isVeri = isSafe([Rend.inf(1:2);Rend.sup(3:4)]);
    % x = [Rend.inf(1:2);Rend.sup(3:4)]';
    
    Q = {};
    for i = 1:4
        Q_i = zeros(4, 4);
        Q_i(i, i) = 1;
        Q{i} = Q_i;
    end

    Rend_isSafe = [v1 0 0 0; 
                   0 v1 0 0; 
                   0 0 1 0; 
                   0 0 0 1] * Rend;
    Rend_isSafe = quadMap(Rend_isSafe, Q);
    Rend_isSafe = [1 1 0 0; 0 0 1 1] * Rend_isSafe;

    evParams_isSafe = struct();
    evParams_isSafe.bound_approx = false;
    evParams_isSafe.reuse_bounds = true;
    evParams_isSafe.add_approx_error_to_GI = true;
    evParams_isSafe.remove_GI = true;
    evParams_isSafe.num_generators = 10000;
    evParams_isSafe.max_bounds = 1;
    % evParams_isSafe.force_approx_lin_at = 0;

    nn_isSafe = neuralNetwork({ ...
        nnRootLayer()} ...
    );
    I = interval(Rend_isSafe, 'split');
    nn_isSafe.layers{1}.l = max([0;0], I.inf);
    nn_isSafe.layers{1}.u = I.sup;
    nn_isSafe.layers{1}.order = 2;
    Rend_isSafe = nn_isSafe.evaluate(Rend_isSafe, evParams_isSafe);

%     Rend_isSafe = [v1 0; 0 1] * Rend_isSafe;
    Rend_isSafe = Rend_isSafe + [v0;0];
    Rend_isSafe = [1,-1] * Rend_isSafe;
    Rend_isSafe = interval(Rend_isSafe, 'split');

    % Rend_isSafe = reduce(Rend_isSafe, 'girard', 500);
    isVeri = 0 <= Rend_isSafe.inf;

    tVeri = toc;
    disp(['Time to check verification: ', num2str(tVeri)]);

    if isVeri
        res = 'VERIFIED';
    else
        res = 'UNKNOWN';
    end
end

tTotal = tSim+tVio+tComp+tVeri;
disp(['Total Time: ', num2str(tTotal)]);
disp(['Result: ' res]);

% Visualization -----------------------------------------------------------

disp("Plotting..")

figure; hold on; box on;
projDims = 1;
useCORAcolors("CORA:contDynamics")
plotOverTime(R, projDims, 'DisplayName', 'Reachable set')
plotOverTime(R(1).R0, projDims, 'DisplayName', 'Initial set');
plotOverTime(simRes, projDims, 'DisplayName', 'Simulations');
xlabel('Time');
ylabel('x');


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
