function completed = example_nonlinear_reach_01_tank
% example_nonlinear_reach_01_tank - example of nonlinear reachability 
%    analysis with conservative linearization
%
% This example can be found in [1, Sec. 3.4.5] or in [2].
%
% Syntax:
%    completed = example_nonlinear_reach_01_tank()
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
% Written:       18-August-2016
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 400;
params.R0 = zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(6)]);
params.U = zonotope([0,0.005]);


% Reachability Settings ---------------------------------------------------

options.timeStep = 1;
options.taylorTerms = 4;
options.zonotopeOrder = 50;
options.alg = 'lin';
options.tensorOrder = 2;

options.lagrangeRem.simplify = 'simplify';


% System Dynamics ---------------------------------------------------------

tank = nonlinearSys(@tank6Eq);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(tank, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 10;
simRes = simulateRandom(tank, params, simOpt);


% Visualization -----------------------------------------------------------

dims = {[1 2],[3 4],[5 6]};
figure;

for k = 1:length(dims)
    
    subplot(1,3,k); hold on; box on;
    projDim = dims{k};
    
    % plot reachable sets
    useCORAcolors("CORA:contDynamics")
    plot(R,projDim,'DisplayName','Reachable set','Unify',true,'UnifyTotalSets',5);
    
    % plot initial set
    plot(R(1).R0,projDim, 'DisplayName','Initial set');
    
    % plot simulation results      
    plot(simRes,projDim,'DisplayName','Simulations');

    % label plot
    xlabel(['x_{',num2str(projDim(1)),'}']);
    ylabel(['x_{',num2str(projDim(2)),'}']);
    legend()
end


% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
