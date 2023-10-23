function res = example_linear_reach_03_inputTrajectory()
% example_linear_reach_03_inputTrajectory - example for linear reachability 
%    analysis with an input trajectory 
%
% Syntax:
%    res = example_linear_reach_03_inputTrajectory()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Matthias Althoff
% Written:       27-July-2018
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

% load model parameters
load data_Jean-MarcBiannic.mat A B uvec
projDims = length(A);
inputDim = length(B(1,:));

% final time, initial set, and uncertain inputs
params.tFinal = 10;
params.R0 = zonotope([zeros(projDims,1),0.1*eye(projDims,4)]);
params.U = zonotope([zeros(inputDim,1),diag([0.05 1])]);

% input trajectory
params.u = uvec(:,1:end-1);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.01;
options.taylorTerms = 4;
options.zonotopeOrder = 50;


% System Dynamics ---------------------------------------------------------

sys = linearSys('JeanMarcSys',A,B);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(sys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simRes = simulateRandom(sys, params);


% Visualization -----------------------------------------------------------

dims = {[1 2],[3 4]};

for k = 1:length(dims)
    
    figure; hold on; legend()
    useCORAcolors("CORA:contDynamics")
    projDims = dims{k};

    % plot reachable sets 
    plot(R,projDims,'Order',10, 'DisplayName','Reachable set');
    
    % plot initial set
    plot(R.R0,projDims, 'DisplayName','Initial set');
    
    % plot simulation results    
    plot(simRes,projDims, 'DisplayName', 'Simulations');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
