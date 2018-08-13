function example_nonlinear_reach_ARCH18_rendezvousSX4p_zonoGirard()
% example_nonlinear_reach_ARCH18_rendezvousSX4p_zonoGirard - example of 
% nonlinear reachability analysis from the ARCH18 friendly competition 
% (instance of rendevouz example). The guard intersections are calculated
% with Girard's method
%
% Syntax:  
%    example_nonlinear_reach_ARCH18_rendezvousSX4p_zonoGirard
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
% Author:       Niklas Kochdumper
% Written:      31-May-2018
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
options.timeStepLoc{2} = 5e-2;
options.timeStepLoc{3} = 2e-1;

options.advancedLinErrorComp = 0;
options.tensorOrder = 2;
options.errorOrder = 5;
options.reductionInterval = inf;
options.loc=1;
options.reductionTechnique='girard';
options.maxError = 150e5*ones(5,1);
options.oldError = zeros(5,1);

options.zonotopeOrder=10; % zonotope order
options.polytopeOrder=3; % polytope order
options.taylorTerms=5;
options.reductionTechnique = 'girard';
options.isHyperplaneMap=0;
options.guardIntersect = 'zonoGirard';
options.enclosureEnables = [1]; % choose enclosure method(s)
options.filterLength = [5,7];

% specify hybrid automata
HA = rendeszvous_nonlinear_passive_hp(); % automatically converted from SpaceEx

% set input:
% get locations
loc = get(HA,'location');
nrOfLoc = length(loc);
for i=1:nrOfLoc
    options.uLoc{i} = 0; % input for simulation
    options.uLocTrans{i} = options.uLoc{i}; % input center for reachability analysis
    options.Uloc{i} = zonotope(0); % input deviation for reachability analysis
end




% Simulation --------------------------------------------------------------

% simulate hybrid automaton
fractionVertices = 0.5;
N = 20;
simRes = cell(N,1);
counter = 1;
for i = 1:N
    if counter < fractionVertices*N
        options.x0 = randPointExtreme(options.R0);
    else
        options.x0 = randPoint(options.R0); 
    end
    
    HA = simulate(HA,options); 
    simRes{i} = get(HA,'trajectory');
    counter = counter + 1;
end




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
verificationCollision = 1;
verificationVelocity = 1;
verificationLOS = 1;
spacecraft = interval([-0.1;-0.1],[0.1;0.1]);

% feasible velocity region as polytope
C = [0 0 1 0 0;0 0 -1 0 0;0 0 0 1 0;0 0 0 -1 0;0 0 1 1 0;0 0 1 -1 0;0 0 -1 1 0;0 0 -1 -1 0];
d = [3;3;3;3;4.2;4.2;4.2;4.2];

% line-of-sight as polytope
Cl = [-1 0 0 0 0;tan(pi/6) -1 0 0 0;tan(pi/6) 1 0 0 0];
dl = [100;0;0];

% passive mode -> check if the spacecraft was hit
for i = 1:length(Rcont{3})    
    temp = interval(Rcont{3}{i}{1});
    if isIntersecting(temp(1:2),spacecraft)
       verificationCollision = 0;
       break;
    end
end

% randezvous attempt -> check if spacecraft inside line-of-sight
for i = 2:length(Rcont{2})  

    temp = interval(Cl*Rcont{2}{i}{1})-dl;

    if any(supremum(temp) > 0)
       verificationLOS = 0;
       break;
    end
end

% randezvous attempt -> check if velocity inside feasible region
for i = 1:length(Rcont{2})  

    temp = interval(C*Rcont{2}{i}{1})-d;

    if any(supremum(temp) > 0)
       verificationVelocity = 0;
       break;
    end
end

verificationCollision
verificationVelocity
verificationLOS

tVer = toc;
disp(['computation time of verification: ',num2str(tVer)]);




% Visualiztion ------------------------------------------------------------

warning('off','all')
figure 
hold on
box on
options.projectedDimensions = [1 2];
options.plotType = {'b','m','g'};
plot(HA,'reachableSet',options); %plot reachable set
plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
for i = 1:length(simRes)
   for j = 1:length(simRes{i}.x)
       plot(simRes{i}.x{j}(:,options.projectedDimensions(1)), ...
            simRes{i}.x{j}(:,options.projectedDimensions(2)),'k'); 
   end
end
xlabel('x');
ylabel('y');

figure 
hold on
grid on
circle(0,0,3);
% fill([-3,-1.2,1.2,3,3,1.2,-1.2,-3],[-1.2,-3,-3,-1.2,1.2,3,3,1.2],'y','FaceAlpha',0.6,'EdgeColor','none');
options.projectedDimensions = [3 4];
options.plotType = {'b','m','g'};
plot(HA,'reachableSetFilled',options); %plot reachable set
plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
for i = 1:length(simRes)
   for j = 1:length(simRes{i}.x)
       plot(simRes{i}.x{j}(:,options.projectedDimensions(1)), ...
            simRes{i}.x{j}(:,options.projectedDimensions(2)),'k'); 
   end
end
warning('on','all')
xlabel('v_x [m/min]');
ylabel('v_y [m/min]');



% Visited Locations -------------------------------------------------------

traj = get(HA,'trajectory');
locTraj = traj.loc;
setTraj = unique(locTraj)

function h = circle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    fill(xunit,yunit,'y','FaceAlpha',0.6,'EdgeColor','none')



%------------- END OF CODE --------------
