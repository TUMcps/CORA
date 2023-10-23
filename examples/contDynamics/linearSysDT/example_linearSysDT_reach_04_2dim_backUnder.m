function completed = example_linearSysDT_reach_04_2dim_backUnder()
% example_linearSysDT_reach_04_2dim_backUnder - example of discrete-time linear 
%    backward reachability analysis with uncertain inputs
%
% Syntax:
%    example_linearSysDT_reach_04_2dim_backUnder
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       21-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% System Dynamics ---------------------------------------------------------

% taken from (82) in "Reach set computation and control synthesis for 
% discrete-time dynamical systems with disturbances"

% system matrix
A = [0 1;...
    -0.5 1];

% input matrix
B = [1; 1];

% constant input
c = zeros(2,1);

% sampling time
dt = 1;

sys = linearSysDT('sys',A,B,c,dt);


% Parameter ---------------------------------------------------------------

params.tFinal = 10;
params.R0 = zonotope([[1;-5],eye(length(A))]);
params.U = zonotope([1,1]);
params.W = [0; 1]*zonotope([0,0.03]);


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 200;
options.reductionTechnique = 'sum';
options.linAlg = 'backward_minmax';


% Reachability Analysis ---------------------------------------------------

tic
R = reach(sys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

% final reachable set becomes new initial set
params.R0 = R.timePoint.set{end};
simOpt.points = 25;
simOpt.type = 'constrained';
simOpt.R = R;
simRes = simulateRandom(sys, params, simOpt);


% Visualization -----------------------------------------------------------

% plot different projections
dims = {[1 2]};

for k = 1:length(dims)
    
    figure; hold on; box on
    useCORAcolors("CORA:contDynamics")
    projDims = dims{k};

    % plot reachable set
    plot(R,projDims);
    
    % plot initial output set
    plot(R(1).R0);
    
    % plot simulation results
    plot(simRes,projDims,'Marker','.');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
