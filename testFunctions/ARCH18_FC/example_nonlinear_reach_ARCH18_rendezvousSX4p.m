function example_nonlinear_reach_ARCH18_rendezvousSX4p()
% example_linear_reach_ARCH18_rendezvousSX4p_simRun - example of linear 
% reachability analysis from the ARCH18 friendly competition 
% (instance of rendevouz example)
%
% Syntax:  
%    example_linear_reach_ARCH18_rendezvousSX4p_simRun
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

% debugging mode
options.debug = 0;

% initial set
R0 = zonotope([[-900; -400; 0; 0; 0],diag([25;25;0;0;0])]);

% initial set
options.x0=center(R0); % initial state for simulation
options.R0=R0; % initial state for reachability analysis

% other
options.startLoc = 1; % initial location
options.finalLoc = -inf; % no final location
options.tStart=0; % start time
options.tFinal=200; % final time
options.intermediateOrder = 2;
options.originContained = 0;
options.timeStepLoc{1} = 2e-1;
options.timeStepLoc{2} = 5e-2;
options.timeStepLoc{3} = 2e-2;

options.advancedLinErrorComp = 0;
options.tensorOrder = 2;
options.errorOrder = 5;
options.reductionInterval = inf;
options.loc=1;
options.reductionTechnique='girard';
options.maxError = 150e5*ones(5,1);
options.oldError = zeros(5,1);

options.zonotopeOrder=40; % zonotope order
options.polytopeOrder=3; % polytope order
options.taylorTerms=3;
options.reductionTechnique = 'girard';
options.isHyperplaneMap=0;
options.guardIntersect = 'polytope';
options.enclosureEnables = [3, 5]; % choose enclosure method(s)
options.filterLength = [5,7];

% specify hybrid automata
HA = rendeszvous_nonlinear_passive(); % automatically converted from SpaceEx
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

%reachable set computations
warning('off','all')
%profile on
tic
[HA] = reach(HA,options);
tComp = toc;
disp(['computation time for spacecraft rendezvous: ',num2str(tComp)]);
%profile viewer
warning('on','all')

% choose projection and plot------------------------------------------------
warning('off','all')
figure 
hold on
grid on
options.projectedDimensions = [1 2];
options.plotType = {'b','m','g'};
plot(HA,'reachableSet',options); %plot reachable set
plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
plot(HA,'simulation',options); % plot simulation
xlabel('x [m]');
ylabel('y [m]');

figure 
hold on
grid on
options.projectedDimensions = [3 4];
options.plotType = {'b','m','g'};
plot(HA,'reachableSet',options); %plot reachable set
plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
plot(HA,'simulation',options); % plot simulation
warning('on','all')
xlabel('v_x [m]');
ylabel('v_y [m/min]');
%--------------------------------------------------------------------------

% detect visited locations-------------------------------------------------
% get results
traj = get(HA,'trajectory');
locTraj = traj.loc;
setTraj = unique(locTraj)
%--------------------------------------------------------------------------


%------------- END OF CODE --------------
