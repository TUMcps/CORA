function testEnv_SXbball1
% test for simulating a bouncing ball specified by an original SX file
%
% Syntax:  
%    testing_env_SXbball
%
% currently only simuation, no reach


%------------- BEGIN CODE --------------

%create hybrid automaton from xml file-------------------------------------
sha = SX2structHA('bball.xml','ballHA','unusedvalue');
computeFile(sha);

HA = ballHA();
%--------------------------------------------------------------------------

%set options---------------------------------------------------------------
options.startLoc = 1; %initial location
options.finalLoc = 0; %0: no final location
options.tStart = 0; %start time
options.tFinal = 1.7; %final time
options.timeStepLoc{1} = 0.05; %time step size for reachable set computation in location 1
options.taylorTerms = 10;
options.zonotopeOrder = 20;
options.polytopeOrder = 10;
options.errorOrder=2;
options.reductionTechnique = 'girard';
options.isHyperplaneMap = 0;
options.enclosureEnables = 5; %choose enclosure method(s)
options.originContained = 0;
%--------------------------------------------------------------------------

%initialize state
options.x0 = [1; 0]; %initial state for simulation
options.R0 = zonotope([options.x0, diag([0.05, 0.05])]); %initial state for reachability analysis

%set input:
options.uLoc{1} = zeros(0,0); %input for simulation
options.uLocTrans{1} = options.uLoc{1}; %input center for reachability analysis
options.Uloc{1} = zonotope(zeros(0,1)); %input deviation for reachability analysis

%simulate hybrid automatons
HA = simulate(HA,options); 

%compute reachable sets
%[HA] = reach(HA,options);

%choose projection and plot------------------------------------------------
figure 
hold on
options.projectedDimensions = [1 2];
options.plotType = 'b';
%plot(HA,'reachableSet',options); %plot reachable set
%plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
plot(HA,'simulation',options); %plot simulation
axis([0,1.2,-6,4]);
%--------------------------------------------------------------------------

%------------- END OF CODE --------------
