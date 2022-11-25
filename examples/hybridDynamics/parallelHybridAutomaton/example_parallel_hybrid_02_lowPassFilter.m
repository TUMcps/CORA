function res = example_parallel_hybrid_02_lowPassFilter()
% example_parallel_hybrid_02_lowPassFilter - example for reachability of a
%    parallel hybrid automaton. The system consists of two piecewise linear
%    low-pass filters that are connected in series
%
% Syntax:  
%    example_parallel_hybrid_02_lowPassFilter
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean, true if completed

% Author:       Niklas Kochdumper
% Written:      06-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------



% System Dynamics ---------------------------------------------------------

PHA = lowpassFilter();


% Parameter ---------------------------------------------------------------

params.tFinal = 0.4;
params.startLoc = [1;3]; 
params.R0 = zonotope([0;0;0;0],diag([0.01,0.01,0.1,0.1]));  


% Reachability Settings ---------------------------------------------------

% settings for continuous reachability 
options.taylorTerms = 8; 
options.zonotopeOrder = 9; 
options.timeStep = 1e-04; 
 
% settings for hybrid systems
options.enclose = {'box'}; 
options.guardIntersect = 'conZonotope';
options.guardOrder = 3;


% Simulation --------------------------------------------------------------

simOpt.points = 10;        % number of initial points
simOpt.fracVert = 0.5;     % fraction of vertices initial set
simOpt.fracInpVert = 0.5;  % fraction of vertices input set
simOpt.inpChanges = 20;    % changes of input over time horizon  

simRes = simulateRandom(PHA,params,simOpt);


% Reachability Analysis ---------------------------------------------------

R = reach(PHA,params,options);


% Visualization -----------------------------------------------------------

% plot filter 1
figure; box on; hold on
plot(R,[1,2],'FaceColor',[.6 .6 .6],'Filled',true,'EdgeColor','none');
plot(simRes,[1,2]);
xlabel('$x_{1}$','interpreter','latex','FontSize',20);
ylabel('$x_{2}$','interpreter','latex','FontSize',20);
title('Filter 1');

% plot filter 2
figure; box on; hold on
plot(R,[3,4],'FaceColor',[.6 .6 .6],'Filled',true,'EdgeColor','none');
plot(simRes,[3,4]);
xlabel('$x_{1}$','interpreter','latex','FontSize',20);
ylabel('$x_{2}$','interpreter','latex','FontSize',20);
title('Filter 2');

res = 1;

%------------- END OF CODE --------------