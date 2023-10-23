function completed = example_nonlinearSysDT_reach_tank
% example_nonlinearSysDT_reach_tank - discrete time veersion of the example
%    of nonlinear reachability analysis with conservative linearization;
%    the continuous-time version can be found in [1, Sec. 3.4.5] or in [2].
%
% Syntax:
%    completed = example_nonlinearSysDT_reach_tank
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] M. Althoff, “Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, 2010
%    [2] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", CDC 2008

% Authors:       Matthias Althoff
% Written:       25-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 100; %final time
params.R0 = zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(6)]);
params.U = zonotope([0,0.005]);


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 50; %zonotope order
options.tensorOrder = 2;
options.errorOrder = 1;
options.lagrangeRem.simplify = 'simplify';


% System Dynamics ---------------------------------------------------------

% sampling time
dt = 0.5;
fun = @(x,u) tank6EqDT(x,u,dt);

tank = nonlinearSysDT('tankSystem',fun,dt); % initialize tank system


% Reachability Analysis ---------------------------------------------------

tic
R = reach(tank, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 60;
simRes = simulateRandom(tank, params, simOpt);


% Visualization -----------------------------------------------------------

dims = {[1 2],[3 4],[5 6]};

for k = 1:length(dims)
    
    figure; hold on; box on;
    projDim = dims{k};
    useCORAcolors("CORA:contDynamics")
    
    % plot reachable sets
    plot(R,projDim);
    
    % plot initial set
    plot(R.R0,projDim);
    
    % plot simulation results     
    plot(simRes,projDim,'Marker','.');

    % label plot
    xlabel(['x_{',num2str(projDim(1)),'}']);
    ylabel(['x_{',num2str(projDim(2)),'}']);
end


% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
