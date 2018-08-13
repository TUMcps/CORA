function example_nonlinear_reach_ARCH18_vanDerPol_zonoGirard()
% example_nonlinear_reach_ARCH18_vanDerPol_zonoGirard - example of
% nonlinear reachability analysis. The guard intersections for the
% pseudoinvariant are calculated with Girards method
%
% Syntax:  
%    example_nonlinear_reach_ARCH18_vanDerPol_zonoGirard()
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
% Written:      02-April-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


%set options---------------------------------------------------------------
options.tStart=0; %start time
%options.tFinal=2; %final time
options.tFinal=7; %final time
options.x0=[1.4; 2.3]; %initial state for simulation

Z0{1}=zonotope([1.4 0.15 0; 2.3 0 0.05]); %initial state for reachability analysis
%Z0{1}=zonotope([1.4 0.6 0; 2.3 0 0.05]); %initial state for reachability analysis
options.R0 = Z0{1};
%options.R0 = zonotopeBundle(Z0);

options.startLoc = 1; %initial location
options.finalLoc = 3; %2: gears are meshed
options.timeStepLoc{1}=0.01; %time step size for reachable set computation
options.timeStepLoc{2}=0.01; %time step size for reachable set computation
options.taylorTerms=3; %number of taylor terms for reachable sets
options.zonotopeOrder=20; %zonotope order
options.polytopeOrder=4; %polytope order

options.advancedLinErrorComp = 0;
options.tensorOrder=1;
options.reductionTechnique='girard';
options.isHyperplaneMap=0;
options.guardIntersect = 'zonoGirard';
options.enclosureEnables = [1];

options.filterLength = [10, 5];
%options.maxError=0.05*[1; 1];
options.maxError=0.5*[1; 1];
options.reductionInterval=inf;
%options.reductionInterval=100;
options.verbose = 1;

%uncertain inputs
options.uTrans = 0;
options.U = zonotope(0); %input for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
vanderPol=nonlinearSys(2,1,@vanderPolEq,options); %initialize van-der-Pol oscillator
%--------------------------------------------------------------------------


%specify hybrid automaton--------------------------------------------------
%define large and small distance
dist = 1e3;
smallDist = 1e3*eps;

%first guard set
guard1 = constrainedHyperplane(halfspace([1, 0],1.5),[0 1;0 -1],[0;2]);
%resets
reset1.A = eye(2); 
reset1.b = zeros(2,1);

%first transition
trans{1} = transition(guard1,reset1,2,'a','b'); %--> next loc: 1; 'a', 'b' are dummies

%invariant
inv = interval([1.5-smallDist; -dist],[dist;dist]);

%specify location
loc{1} = location('loc1',1,inv,trans,vanderPol); 
loc{2} = location('loc2',2,ones(2,1)*interval(-dist,dist),[],vanderPol); 
%specify hybrid automata
HA = hybridAutomaton(loc); % for "geometric intersection"
%--------------------------------------------------------------------------

%set input:
for i=1:2
    options.uLoc{i} = 0; %input for simulation
    options.uLocTrans{i} = options.uLoc{i}; %input center for reachability analysis
    options.Uloc{i} = zonotope(0); %input deviation for reachability analysis
end

%obtain random simulation results
for i=1:10
    %set initial state, input
    if i <=4 
        options.x0 = randPointExtreme(options.R0); %initial state for simulation
    else
        options.x0 = randPoint(options.R0); %initial state for simulation
    end 
    
    %simulate hybrid automaton
    HAsim = simulate(HA,options); 
    simRes{i} = get(HAsim,'trajectory');
end

HA = simulate(HA,options); 

%reachable set computations
%profile on
tic
[HA] = reach(HA,options);
tComp = toc;
disp(['computation time for van der Pol: ',num2str(tComp)]);
%profile off
%profile viewer

%choose projection and plot------------------------------------------------
figure 
hold on
Rcont = get(HA,'continuousReachableSet');
%plot continuous reachable set
redOrder = 4;
%1st time interval
Rplot = Rcont.OT;
%Rplot_noHP = Rcont_noHP.OT;
for iMode = 1:length(Rplot)
    for iSet = 1:length(Rplot{iMode})
        Rtmp = reduce(project(Rplot{iMode}{iSet}{1},[1 2]),'girard',redOrder);
        plotFilled(Rtmp,[1 2],[.6 .6 .6],'EdgeColor','none');            
    end
end

plotFilled(options.R0,[1 2],'w','EdgeColor','k'); %plot initial set

%plot simulation results
for n=1:length(simRes)
    for j=1:length(simRes{n}.x)
        %plot results
        plot(simRes{n}.x{j}(:,1), simRes{n}.x{j}(:,2),'k');
    end
end

xlabel('x');
ylabel('y');
box on
%--------------------------------------------------------------------------




%------------- END OF CODE --------------