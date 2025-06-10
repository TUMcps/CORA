function [completed,res,tTotal] = example_neuralNet_reach_10_spacecraftDocking
% example_neuralNet_reach_10_spacecraftDocking - example of reachability analysis
%    for an neural network controlled system
%
% Syntax:
%    completed = example_neuralNet_reach_10_spacecraftDocking()
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
% Written:       20-June-2022
% Last update:   30-March-2023 (TL, ARCH'23 revisions)
%                06-May-2025 (TL, ARCH'25 revisions, nn not verifiably robust..)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp("BENCHMARK: Spacecraft Docking")

% Parameter ---------------------------------------------------------------

m = 12;
n = 0.001027;

params.tFinal = 40;
params.R0 = polyZonotope(interval( ...
    [70; 70; -0.28; -0.28], ...
    [106; 106; 0.28; 0.28] ...
));
params.R0 = params.R0;

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.1;
options.alg = 'lin';
options.tensorOrder = 2;
options.taylorTerms = 4;
options.zonotopeOrder = 20;

% Options for NN evaluation -----------------------------------------------

options.nn = struct();
options.nn.poly_method = "regression";

% System Dynamics ---------------------------------------------------------

% open-loop system
% f(x, u) = A*x + B*u
A = [0 0 1 0; 
     0 0 0 1; 
     3*n^2 0 0 2*n; 
     0 0 -2*n 0];
B = [0 0; 
     0 0; 
     1/m 0; 
     0 1/m;];
f = @(x,u) A*x(1:4,:) + B*u(1:2,:);
sys = nonlinearSys(f);

% load neural network controller
% [4, 256, 256, 2]
nn = neuralNetwork.readONNXNetwork('controller_spacecraftDocking.onnx');

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.1);

% Specification -----------------------------------------------------------

isSafe = @(x) (sqrt(x(3)^2 + x(4)^2)) <= (0.2 + 2*n*sqrt(x(1)^2 + x(2)^2));

% Simulation --------------------------------------------------------------

timerVal = tic;
simRes = simulateRandom(sys, params);
tSim = toc(timerVal);
disp(['Time to compute random simulations: ', num2str(tSim)]);

% Check Violation --------------------------------------------------------

timerVal = tic;
isVio = aux_checkSimulations(simRes,isSafe);
tVio = toc(timerVal);
disp(['Time to check violation in simulations: ', num2str(tVio)]);


if isVio 
    res = 'VIOLATED';
    R = [];
    tComp = 0;
    tVeri = 0;
else
    % Reachability Analysis -----------------------------------------------

    timerVal = tic;
    R = reach(sys, params, options);
    tComp = toc(timerVal);
    disp(['Time to compute reachable set: ', num2str(tComp)]);

    % Verification --------------------------------------------------------

    timerVal = tic;
    isVeri = aux_checkReachSet(R,n,m);
    tVeri = toc(timerVal);
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
projDims = 1; % [3,4] look more problematic even after single step
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


% Auxiliary functions -----------------------------------------------------

function isVio = aux_checkSimulations(simRes,isSafe)
    % check specifications in simulations
    
    % iterate through all simulations
    for i = 1:numel(simRes)
        x_i = simRes(i).x;
        for j = 1:numel(x_i)
            x_ij = x_i{j};
            for k = 1:size(x_ij,1)
                % check each point
                isVio = ~isSafe(x_ij(k, :)');
                if isVio
                    return;
                end
            end
        end
    end

    % no violating point found
    isVio = false;
end

function isVeri = aux_checkReachSet(R,n,m)
    % check specification in reachable set

    % iterate through all timeInterval solutions
    for i=1:numel(R)
        R_i = R(i);
        for j=1:numel(R_i.timeInterval.set)
            % get current reachable set
            R_ij = polyZonotope(R_i.timeInterval.set{j});
            time_ij = R_i.timeInterval.time{j};

            % apply order reduction and restructure GI
            R_ij = reduce(R_ij,'girard',4);
            R_ij = polyZonotope( ...
                R_ij.c, ...
                [R_ij.G R_ij.GI], ...
                [],blkdiag(R_ij.E,eye(size(R_ij.GI,2))));

            % compute x^2
            Q = cell(1,4);
            for q = 1:4
                Q_i = zeros(4, 4);
                Q_i(q, q) = 1;
                Q{q} = Q_i;
            end
            R_ij = quadMap(R_ij, Q);

            % sum respective terms
            R_ij = [1 1 0 0; 0 0 1 1] * R_ij;

            % compute 2*n before computing root to reduce domain
            R_ij = [(2*n)^2 0;0 1] * R_ij;

            % apply order reduction
            % R_ij = reduce(R_ij,'girard',4);
        
            % compute square root
            nn_sqrt = neuralNetwork({ ...
                nnRootLayer()} ...
            );
            R_ij_int = interval(R_ij);
            nn_sqrt.layers{1}.l = max(R_ij_int.inf, [0;0]);
            nn_sqrt.layers{1}.u = R_ij_int.sup;
            nn_sqrt.layers{1}.order = 1;
            optionsSqrt = struct;
            optionsSqrt.nn.remove_GI = false;
            optionsSqrt.nn.reuse_bounds = true;
            optionsSqrt.nn.num_generators = 10000;
            optionsSqrt.nn.poly_method = 'singh';
            R_ij = nn_sqrt.evaluate(R_ij,optionsSqrt);
        
            % finish each side of inequality
            R_ij = R_ij + [0.2;0];

            % determine safety
            R_ij = [1,-1] * R_ij;
            R_ij_int = interval(R_ij);
            isVeri = 0 <= R_ij_int.inf;
            if ~isVeri
                % try again with splitting
                R_ij_int = interval(R_ij, 'split');
                isVeri = 0 <= R_ij_int.inf;
            end

            if ~isVeri
                fprintf('Unable to verify t \\in [%.2f,%.2f].\n',time_ij.inf,time_ij.sup);
                return
            end
        end
    end

    % no violation found
    isVeri = true;   
end

% ------------------------------ END OF CODE ------------------------------
