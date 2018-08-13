function example_linear_reach_ARCH18_rendezvousSX4np_simRun()
% example_linear_reach_ARCH18_rendezvousSX4np_simRun - example of linear 
% reachability analysis from the ARCH18 friendly competition 
% (instance of rendevouz example)
%
% Syntax:  
%    example_linear_reach_ARCH18_rendezvousSX4np_simRun
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      20-April-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% x_1 = v_x
% x_2 = v_y
% x_3 = s_x 
% x_4 = s_y
% x_5 = t


% initial set
R0 = zonotope([[-900; -400; 0; 0; 0],[25;25;0;0; 0]]);

% initial set
options.x0=center(R0); % initial state for simulation
options.R0=R0; % initial state for reachability analysis

% other
options.startLoc = 1; % initial location
options.finalLoc = -inf; % no final location
options.tStart=0; % start time
options.tFinal=200; % final time
% options.tFinal=0.5; % final time
% options.tFinal=inf; % final time
options.intermediateOrder = 2;
options.originContained = 0;
options.timeStepLoc{1} = 1e-3;
options.timeStepLoc{2} = 1e-3;

options.zonotopeOrder=40; % zonotope order
options.polytopeOrder=3; % polytope order
options.taylorTerms=3;
options.reductionTechnique = 'girard';
options.isHyperplaneMap=0;
options.enclosureEnables = [1, 2, 3, 5]; % choose enclosure method(s)
options.filterLength = [5,7];

% specify hybrid automata
HA = rendezvousSX4np(); % automatically converted from SpaceEx
%--------------------------------------------------------------------------

% set input:
% get locations
loc = get(HA,'location');
nrOfLoc = length(loc);
for i=1:nrOfLoc
    options.uLoc{i} = 0; % input for simulation
    options.uLocTrans{i} = options.uLoc{i}; % input center for reachability analysis
    options.Uloc{i} = zonotope(0); % input deviation for reachability analysis
end

% simulate hybrid automaton
HA = simulate(HA,options); 


% choose projection and plot------------------------------------------------
figure 
hold on
options.projectedDimensions = [1 2];
plot(HA,'simulation',options); % plot simulation

figure 
hold on
options.projectedDimensions = [3 4];
plot(HA,'simulation',options); % plot simulation
%--------------------------------------------------------------------------

% detect visited locations-------------------------------------------------
% get results
traj = get(HA,'trajectory');
locTraj = traj.loc;
setTraj = unique(locTraj)
%--------------------------------------------------------------------------


%------------- END OF CODE --------------
