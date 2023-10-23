function completed = example_linearSysDT_reach_01_5dim()
% example_linearSysDT_reach_01_5dim - example of discrete-time linear 
%    reachability analysis with uncertain inputs
%
% Syntax:
%    example_linearSysDT_reach_01_5dim
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       20-March-2020
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% System Dynamics ---------------------------------------------------------

% system matrix
A = [0.95 -0.15 0 0 0;...
    0.15 0.95 0 0 0; ...
    0 0 0.9 0.05 0; ...
    0 0 -0.05 0.9 0; ...
    0 0 0 0 0.92];

% input matrix
B = 1;

% constant input
c = zeros(5,1);

% output matrix
C = [2 0 0 0 0; 0 1 0 0 0];

% sampling time
dt = 0.04;

fiveDimSys = linearSysDT('fiveDimSys',A,B,c,C,dt);


% Parameter ---------------------------------------------------------------

params.tFinal = 5;
params.R0 = zonotope([ones(length(A),1),0.1*eye(length(A))]);
params.U = zonotope([zeros(5,1),0.02*diag([0.1, 0.3, 0.1, 0.3, 0.3])]);


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 200;


% Reachability Analysis ---------------------------------------------------

tic
R = reach(fiveDimSys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 25;
simRes = simulateRandom(fiveDimSys, params, simOpt);


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
    plot(R.timePoint.set{1},projDims);
    
    % plot simulation results
    plot(simRes,projDims,'Traj','y','Marker','.');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
