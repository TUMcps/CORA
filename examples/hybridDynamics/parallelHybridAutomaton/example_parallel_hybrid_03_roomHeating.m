function res = example_parallel_hybrid_03_roomHeating()
% example_parallel_hybrid_03_roomHeating - example for reachability of a
%    parallel hybrid automaton considering the room heating benchmark
%    described in Sec. 2.3 in [1] with two rooms
%
% Syntax:  
%    example_parallel_hybrid_03_roomHeating
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean, true if completed
%
% References:
%   [1] A. Fehnker and F. Ivancic. "Benchmarks for Hybrid Systems 
%       Verification", HSCC 2004

% Author:       Niklas Kochdumper
% Written:      26-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% System Dynamics ---------------------------------------------------------

PHA = roomHeatingParallel();


% Parameter ---------------------------------------------------------------

params.tFinal = 12;
params.startLoc = [1;1]; 
params.R0 = zonotope([20.5;20.5],diag([0.1,0.1]));  
params.U = zonotope(4,0.01);


% Reachability Settings ---------------------------------------------------

% settings for continuous reachability 
options.taylorTerms = 5; 
options.zonotopeOrder = 10; 
options.timeStep = 0.005; 
 
% settings for hybrid systems
options.enclose = {'box','pca'}; 
options.guardIntersect = 'zonoGirard';


% Simulation --------------------------------------------------------------

simOpt.points = 10;        % number of initial points
simOpt.fracVert = 0.5;     % fraction of vertices initial set
simOpt.fracInpVert = 0.5;  % fraction of vertices input set
simOpt.inpChanges = 20;    % changes of input over time horizon  

simRes = simulateRandom(PHA,params,simOpt);


% Reachability Analysis ---------------------------------------------------

R = reach(PHA,params,options);


% Visualization -----------------------------------------------------------

% temperature room 1 over time
figure; box on; hold on
plotOverTime(R,1,'FaceColor',[.6 .6 .6],'EdgeColor','none');
plotOverTime(simRes,1);
xlabel('Time')
ylabel('Temperature');
title('Room 1');
xlim([0,params.tFinal]);

% temperature room 2 over time
figure; box on; hold on
plotOverTime(R,2,'FaceColor',[.6 .6 .6],'EdgeColor','none');
plotOverTime(simRes,2);
xlabel('Time')
ylabel('Temperature');
title('Room 2');
xlim([0,params.tFinal]);

res = 1;

%------------- END OF CODE --------------