function example_linear_reach_ARCH18_rendezvousSX4np_zonoGirard()
% example_linear_reach_ARCH18_rendezvousSX4np_zonoGirard - example of linear 
% reachability analysis from the ARCH18 friendly competition 
% (instance of rendevouz example). The guard intersections are cacluated
% with the Girard's method. 
%
% Syntax:  
%    example_linear_reach_ARCH18_rendezvousSX4np_zonoGirard
%
% Inputs:
%    no
%
% Outputs:
%    no
%
% References: 
%   [1] N. Chan et al. "Verifying safety of an autonomous spacecraft 
%       rendezvous mission (Benchmark proposal)"  

% Author:       Niklas Kochdumper
% Written:      23-May-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% x_1 = v_x
% x_2 = v_y
% x_3 = s_x 
% x_4 = s_y
% x_5 = t



% Options -----------------------------------------------------------------

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
options.timeStepLoc{2} = 2e-2;

options.zonotopeOrder=10; % zonotope order
options.polytopeOrder=3; % polytope order
options.taylorTerms=3;
options.reductionTechnique = 'girard';
options.isHyperplaneMap=0;
options.enclosureEnables = [1]; % choose enclosure method(s)
options.filterLength = [5,7];
options.guardIntersect = 'zonoGirard';

% specify hybrid automata
HA = rendezvousSX4np_hs(); % automatically converted from SpaceEx



% Simulation --------------------------------------------------------------

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



% Reachability Analysis ---------------------------------------------------

%reachable set computations
warning('off','all')
%profile on
tic
[HA] = reach(HA,options);
tComp = toc;
disp(['computation time for spacecraft rendezvous: ',num2str(tComp)]);
%profile viewer
warning('on','all')



% Verification ------------------------------------------------------------

tic
Rcont = get(HA,'continuousReachableSet');
Rcont = Rcont.OT;
verificationVelocity = 1;
verificationLOS = 1;

% feasible velocity region as polytope
C = [0 0 1 0 0;0 0 -1 0 0;0 0 0 1 0;0 0 0 -1 0;0 0 1 1 0;0 0 1 -1 0;0 0 -1 1 0;0 0 -1 -1 0];
d = [3;3;3;3;4.2;4.2;4.2;4.2];

% line-of-sight as polytope
Cl = [-1 0 0 0 0;tan(pi/6) -1 0 0 0;tan(pi/6) 1 0 0 0];
dl = [100;0;0];


% randezvous attempt -> check if spacecraft inside line-of-sight
for i = 2:length(Rcont{2})  

    temp = interval(Cl*Rcont{2}{i})-dl;

    if any(supremum(temp) > 0)
       verificationLOS = 0;
       break;
    end
end

% randezvous attempt -> check if velocity inside feasible region
for i = 1:length(Rcont{2})  

    temp = interval(C*Rcont{2}{i})-d;

    if any(supremum(temp) > 0)
       verificationVelocity = 0;
       break;
    end
end

verificationVelocity
verificationLOS

tVer = toc;
disp(['computation time of verification: ',num2str(tVer)]);




% Visualization -----------------------------------------------------------

% Figure 4 (a), (b) and (c) in reference paper [1]
warning('off','all')
figure 
hold on
box on
fill([-100,0,-100],[-60,0,60],'y','FaceAlpha',0.6,'EdgeColor','none');
options.projectedDimensions = [1 2];
options.plotType = {'b','m'};
plot(HA,'reachableSet',options); %plot reachable set
plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
plot(HA,'simulation',options); % plot simulation
xlabel('s_x');
ylabel('s_y');

% Figure 4 (d) in reference paper[1]
figure 
hold on
grid on
fill([-3,-1.2,1.2,3,3,1.2,-1.2,-3],[-1.2,-3,-3,-1.2,1.2,3,3,1.2],'y','FaceAlpha',0.6,'EdgeColor','none');
options.projectedDimensions = [3 4];
options.plotType = {'b','m'};
plot(HA,'reachableSet',options); %plot reachable set
plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
plot(HA,'simulation',options); % plot simulation
warning('on','all')
xlabel('v_x [m/min]');
ylabel('v_y [m/min]');



% Visited Locations -------------------------------------------------------

% get results
traj = get(HA,'trajectory');
locTraj = traj.loc;
setTraj = unique(locTraj)


%------------- END OF CODE --------------
