function example_nonlinHybrid_reach_ARCH18_pancreas2()
% example_nonlinHybrid_reach_ARCH18_pancreas - example of hybrid 
% nonlinear reachability analysis from the ARCH18 friendly competition.
% Girard's method is used to calculate the intersections with the guards
%
% Syntax:  
%    example_nonlinHybrid_reach_ARCH18_pancreas
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% References:
%   [1] X. Chen et al. "Formal Verification of a Multi-Basal Insulin 
%       Infusion Control Model"

% Author:       Niklas Kochdumper
% Written:      22-May-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% x_1 = X
% x_2 = I_sc1
% x_3 = I_sc2
% x_4 = G_t
% x_5 = G_p
% x_6 = I_l
% x_7 = I_p
% x_8 = I_1
% x_9 = I_d
% x_10 = G_s
% x_11 = t


% Options -----------------------------------------------------------------

% debugging mode
options.debug = 0;

% initial set
int = interval([0;72.429;141.149;162.449;229.824;3.19;5.49;100.249;100.249;120;0;50], ...
               [0;72.431;141.151;162.451;268.128;3.21;5.51;100.251;100.251;160;0;90]);
           
R0 = zonotope(int);

options.x0=center(R0); % initial state for simulation
options.R0=R0; % initial state for reachability analysis

% additional options
options.startLoc = 3; % initial location
options.finalLoc = -inf; % no final location
options.tStart=0; % start time
options.tFinal=719.9; % final time
options.tFinal = 42;
options.intermediateOrder = 10;
options.originContained = 0;
options.timeStepLoc = repmat({0.5},[30,1]);
% options.timeStepLoc{9} = 0.05;

options.zonotopeOrder=10; % zonotope order
options.polytopeOrder=3; % polytope order
options.reductionTechnique = 'girard';
options.isHyperplaneMap=0;
options.guardIntersect = 'polytope';     % use constrained zonotopes to calculate guard intersections
options.intersectInvariant = 0; 
options.enclosureEnables = [3,5]; % choose enclosure method(s)
options.filterLength = [5,7];

options.taylorTerms=20;             % number of taylor terms for reachable sets
options.reachabilitySteps=100;      % number of reachability steps 
options.timeStep=options.tFinal/options.reachabilitySteps; % length of each time step in reachability analysis

options.advancedLinErrorComp = 0;
options.tensorOrder = 1;
options.errorOrder = 10;
options.reductionInterval = inf;
options.loc=3;
options.reductionTechnique='girard';
options.maxError = 150e5*ones(12,1);
options.oldError = zeros(12,1);

% specify hybrid automata
HA = pancreas2(); % automatically converted from SpaceEx

% set input:
loc = get(HA,'location');
nrOfLoc = length(loc);
for i=1:nrOfLoc
    options.uLoc{i} = 0; % input for simulation
    options.uLocTrans{i} = options.uLoc{i}; % input center for reachability analysis
    options.Uloc{i} = zonotope(0); % input deviation for reachability analysis
end




% Simulation --------------------------------------------------------------

for i=1:10
    %set initial state, input
    if i < 5
        options.x0 = randPointExtreme(options.R0);
    else
        options.x0 = randPoint(options.R0);
    end 
    
    %simulate hybrid automaton
    HAsim = simulate(HA,options); 
    simRes{i} = get(HAsim,'trajectory');
end

% simulate hybrid automaton
HA = simulate(HA,options);



% Reachability Analysis ---------------------------------------------------

warning('off','all')
% profile on
tic
[HA] = reach(HA,options);
tComp = toc;
disp(['computation time for artificial pancreas: ',num2str(tComp)]);
% profile viewer
warning('on','all')



% Visualization -----------------------------------------------------------

% warning('off','all')
% figure 
% hold on
% grid on
% options.plotType = {'b','m','g','r','y','b','m','g','r','y','b','m','g','r','y','b','m','g','r','y','b','m','g','r','y','b','m','g','r','y'};
% options.projectedDimensions = [1 2];
% plot(HA,'reachableSet',options); %plot reachable set
% plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
% plot(HA,'simulation',options); % plot simulation
% xlabel('X');
% ylabel('Isc_1');
% 
% figure 
% hold on
% grid on
% options.projectedDimensions = [3 4];
% plot(HA,'reachableSet',options); %plot reachable set
% plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
% plot(HA,'simulation',options); % plot simulation
% warning('on','all')
% xlabel('Isc_2');
% ylabel('G_t');
% 
% figure 
% hold on
% grid on
% options.projectedDimensions = [5 6];
% plot(HA,'reachableSet',options); %plot reachable set
% plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
% plot(HA,'simulation',options); % plot simulation
% warning('on','all')
% xlabel('G_p');
% ylabel('I_l');
% 
% figure 
% hold on
% grid on
% options.projectedDimensions = [7 8];
% plot(HA,'reachableSet',options); %plot reachable set
% plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
% plot(HA,'simulation',options); % plot simulation
% warning('on','all')
% xlabel('I_p');
% ylabel('I_1');
% 
% figure 
% hold on
% grid on
% options.projectedDimensions = [9 10];
% plot(HA,'reachableSet',options); %plot reachable set
% plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
% plot(HA,'simulation',options); % plot simulation
% warning('on','all')
% xlabel('I_d');
% ylabel('G_s');

figure 
hold on
box on
options.plotType = {'b','m','g','r','y','b','m','g','r','y','b','m','g','r','y','b','m','g','r','y','b','m','g','r','y','b','m','g','r','y'};
options.projectedDimensions = [11 10];
plot(HA,'reachableSet',options); %plot reachable set
plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
for n=1:length(simRes)
    for j=1:length(simRes{n}.x)
        plot(simRes{n}.x{j}(:,options.projectedDimensions(1)), simRes{n}.x{j}(:,options.projectedDimensions(2)),'k');
    end
end
warning('on','all')
xlabel('t');
ylabel('G_s');
xlim([0,720]);



% Visited Locations -------------------------------------------------------

traj = get(HA,'trajectory');
locTraj = traj.loc;
setTraj = unique(locTraj)


%------------- END OF CODE --------------