function completed = example_nonlinearDA_reach_01_powerSystem_3bus()
% example_nonlinearDA_reach_01_powerSystem_3bus - example of 
%    nonlinear-differential-algebraic reachability analysis, 3-bus power system
%
% Syntax:  
%    example_nonlinearDA_reach_01_powerSystem_3bus()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% References:
%    -

% Author:       Matthias Althoff
% Written:      18-August-2016
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---


%------------- BEGIN CODE --------------

% Parameter ---------------------------------------------------------------

nrOfConstr = 6;
params.tFinal = 5; 

x0 = [380; 0.7];
params.y0guess = [ones(0.5*nrOfConstr, 1); zeros(0.5*nrOfConstr, 1)];
params.R0 = zonotope([x0,diag([0.1, 0.01])]);

params.U = zonotope([[1; 0.4],diag([0, 0.04])]);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.05;
options.taylorTerms = 6;
options.zonotopeOrder = 200;
options.tensorOrder = 2;

options.maxError_x = [0.5; 0];
options.maxError_y = 0.005*[1; 1; 1; 1; 1; 1];


% System Dynamics ---------------------------------------------------------

powerDyn = nonlinDASys(@bus3Dyn,@bus3Con);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(powerDyn, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 60;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;

simRes = simulateRandom(powerDyn, params, simOpt);


% Visualization -----------------------------------------------------------

dim = [1 2];
    
figure; hold on; box on;

% plot reachable sets
plot(R,dim,'FaceColor',[.7 .7 .7],'EdgeColor','none');

% plot initial set
plot(params.R0,dim,'w','Filled',true,'EdgeColor','k');

% plot simulation results      
plot(simRes,dim);

% label plot
xlabel(['x_{',num2str(dim(1)),'}']);
ylabel(['x_{',num2str(dim(2)),'}']);

% example completed
completed = 1;

%------------- END OF CODE --------------
        