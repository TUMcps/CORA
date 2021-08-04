function completed = example_nonlinear_reach_01_tank
% example_nonlinear_reach_01_tank - example of nonlinear reachability 
%                                   analysis with conservative linearizaton
%
% This example can be found in [1, Sec. 3.4.5] or in [2].
%
% Syntax:  
%    completed = example_nonlinear_reach_01_tank()
%
% Inputs:
%    no
%
% Outputs:
%    completed - boolean 
%
% References:
%    [1] M. Althoff, “Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, 2010
%    [2] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", CDC 2008

% Author:       Matthias Althoff
% Written:      18-August-2016
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---

%------------- BEGIN CODE --------------

% Parameter ---------------------------------------------------------------

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

simOpt.points = 60;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;

simRes = simulateRandom(tank, params, simOpt);


% Visualization -----------------------------------------------------------

dims = {[1 2],[3 4],[5 6]};

for k = 1:length(dims)
    
    figure; hold on; box on;
    projDim = dims{k};
    
    % plot reachable sets
    plot(R,projDim,'FaceColor',[.8 .8 .8],'EdgeColor','none');
    
    % plot initial set
    plot(params.R0,projDim,'w','Filled',true,'EdgeColor','k');
    
    % plot simulation results     
    plot(simRes,projDim,'b');

    % label plot
    xlabel(['x_{',num2str(projDim(1)),'}']);
    ylabel(['x_{',num2str(projDim(2)),'}']);
end


% example completed
completed = 1;

%------------- END OF CODE --------------
