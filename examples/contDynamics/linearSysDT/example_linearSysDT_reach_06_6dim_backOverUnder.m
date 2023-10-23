function completed = example_linearSysDT_reach_06_6dim_backOverUnder()
% example_linearSysDT_reach_06_6dim_backOverUnder - example of 
%    discrete-time linear backward reachability analysis with uncertain 
%    inputs taken from Sec. VI.A in [1].
%    This example considers the longitudinal dynamics. 
%
% Syntax:
%    example_linearSysDT_reach_06_6dim_backOverUnder
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
    [0.9911, -0.04858, -0.01709, -0.4883, 0, 0; ...
    0.0005870, 0.9968, 0.5168, -0.0001398, 0, 0; ...
    0.0002070, -0.001123, 0.9936, -5.092e-5, 0, 0; ...
    1.907, -1.032, 0.01832, 1, 0, 0; ...
    -0.04601, 0.001125, 0.0002638, 0.01130, 1, 0; ...
    -5.095e-5, -0.1874, -0.01185, 4.004, 0, 1];

% input matrix 
B = [...
    1.504, 7.349e-5; ...
    -0.04645, -3.421e-6; ...
    -0.009812, -1.488e-6; ...
    -9.080e-5, -1.371e-8; ...
    -0.03479, -1.700e-6; ...
    0.004171, 2.913e-7];

% constant input
c = zeros(length(A),1);

% instantiate linear discrete time system
sys = linearSysDT('sys',A,B,c,dt);


% Parameter ---------------------------------------------------------------

%params.tFinal = 80;
params.tFinal = 4;
params.R0 = zonotope([[50;5;0;0;400;325],diag([10, 5, 0.1, pi, 400, 65])]);
params.U = zonotope([[0.131; 5e3], diag([0.393, 5e3])]);
params.W = zonotope([zeros(6,1),diag([0.3025, 0.4025, 0.01213, 0.006750, 1.373, 1.331])]);


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 3;


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
    plot(R_over(1).R0,projDims);
    
    % plot simulation results
    plot(simRes,projDims,'Marker','.');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
