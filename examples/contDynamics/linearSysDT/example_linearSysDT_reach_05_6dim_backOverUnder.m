function completed = example_linearSysDT_reach_05_6dim_backOverUnder()
% example_linearSysDT_reach_05_6dim_backOverUnder - example of 
%    discrete-time linear backward reachability analysis with uncertain 
%    inputs taken from Sec. VI.A in [1].
%    This example considers the lateral dynamics. 
%
% Syntax:
%    example_linearSysDT_reach_05_6dim_backOverUnder
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] L. Yang and N. Ozay, "Scalable Zonotopic Under-Approximation of 
%        Backward Reachable Sets for Uncertain Linear Systems," in IEEE 
%        Control Systems Letters, vol. 6, pp. 1555-1560, 2022.

% Authors:       Matthias Althoff
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% System Dynamics ---------------------------------------------------------

% sampling time
dt = 1;

% system matrix
A = ...
    [1.004, 0.1408, 0.3095, -0.3112, 0, 0; ...
    0.03015, 1.177, 0.6016, -0.6029, 0, 0; ...
    -0.02448, -0.1877, 0.3803, 0.5642, 0, 0; ...
    -0.01057, -0.09588, -0.3343, 1.277, 0, 0; ...
    0.0003943, 0.0095901, -0.005341, -0.007447, 1, 0; ...
    -0.2579, -23.32, -51.03, 61.35, -37.86, 1];

% input matrix 
B = [...
    -0.1189, 0.007812; ...
    -0.1217, 0.2643; ...
    0.01773, -0.2219; ...
    -0.02882, -0.09982; ...
    -0.0005607, 0.002437; ...
    0.1120, -0.5785];

% constant input
c = zeros(length(A),1);

% instantiate linear discrete time system
sys = linearSysDT('sys',A,B,c,dt);


% Parameter ---------------------------------------------------------------

params.tFinal = 40;
params.R0 = zonotope([zeros(6,1),diag([1, 1, 1, pi/5, pi/5, 2])]);
params.U = zonotope([zeros(2,1), pi*eye(2)]);
params.W = zonotope([zeros(6,1),diag([0.037, 0.00166, 0.0078, 0.00124, 0.00107, 0.07229])]);


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 6;


% Reachability Analysis (over-approximative) ------------------------------

options.linAlg = 'backward_maxmin_coarse';

tic
R_over = reach(sys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

% Reachability Analysis (under-approximative) -----------------------------

options.reductionTechnique = 'sum';
%options.linAlg = 'backward_minmax_RKweighted';
options.linAlg = 'backward_minmax';

tic
R_under = reach(sys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

% final reachable set becomes new initial set
params.R0 = R_under.timePoint.set{end};
simOpt.points = 25;
simOpt.type = 'constrained';
simOpt.R = R_under; % requires under-approximative reachable sets
simRes = simulateRandom(sys, params, simOpt);


% Visualization -----------------------------------------------------------

% plot different projections
dims = {[1 2]};

for k = 1:length(dims)
    
    figure; hold on; box on
    useCORAcolors("CORA:contDynamics")
    projDims = dims{k};

    % plot reachable set
    plot(R_over,projDims);
    
    % plot initial output set
    plot(R_over.R0,projDims,'FaceColor','w','EdgeColor','k');
    
    % plot simulation results
    plot(simRes,projDims,'Marker','.');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
