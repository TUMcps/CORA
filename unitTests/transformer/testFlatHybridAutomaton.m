function HA = testFlatHybridAutomaton()

% Constant
a = 7.9*10^-5;
% Text = 273.15; % 0K
Text = 0;
P1 = 210*10^3;
P2 = 160*10^3;
C = 7.6*10^7;

temperatureInterval = [Text+10,Text+30;Text+10,Text+30];

% Location 1A : Heater (on,on)
A = [-2*a,a;a,-2*a];
B = [a*Text + P1 / C,0;0,a*Text + P2 / C];

linSys = linearSys('lin',A,B);

inv = intervalhull(temperatureInterval);

reset.A = eye(2);
reset.b = 0;

guard = halfspace([-1,0]',-(Text+22));

trans{1} = transition(guard,reset,4,'a','b'); % to off/on

guard = halfspace([0,-1]',-(Text+18));

trans{2} = transition(guard,reset,2,'a','b'); % to on/off

guard = mptPolytope([-1,0;0,-1],[-(Text+22);-(Text+18)]);

trans{3} = transition(guard,reset,3,'a','b'); % to off/off

loc{1} = location('on/on',1,inv,trans,linSys);

% Location 1B : Heater (on, off)
A = [-2*a,a;a,-2*a];
B = [a*Text + P1 / C,0;0,a*Text];

linSys = linearSys('lin',A,B);

inv = intervalhull(temperatureInterval);

reset.A = eye(2);
reset.b = 0;

guard = halfspace([0,1]',Text+14);

trans{1} = transition(guard,reset,1,'a','b'); % to on/on

guard = halfspace([-1,0]',-(Text+22));

trans{2} = transition(guard,reset,3,'a','b'); % to off/off

guard = mptPolytope([-1,0;0,1],[-(Text+22);Text+14]);

trans{3} = transition(guard,reset,4,'a','b'); % to off/on

loc{2} = location('on/off',4,inv,trans,linSys);

% Location 2B : Heater (off,off)
A = [-2*a,a;a,-2*a];
B = [a*Text,0;0,a*Text];

linSys = linearSys('lin',A,B);

inv = intervalhull(temperatureInterval);

reset.A = eye(2);
reset.b = 0;

guard = halfspace([1,0]',Text+18);

trans{1} = transition(guard,reset,2,'a','b'); % to on/off

guard = halfspace([0,1]',Text+14);

trans{2} = transition(guard,reset,4,'a','b'); % to off/on

guard = mptPolytope([1,0;0,1],[Text+18;Text+14]);

trans{3} = transition(guard,reset,1,'a','b'); % to on/on

loc{3} = location('on/off',3,inv,trans,linSys);

% Location 2A : Heater (off,on)
A = [-2*a,a;a,-2*a];
B = [a*Text,0;0,a*Text + P2 / C];

linSys = linearSys('lin',A,B);

inv = intervalhull(temperatureInterval);

reset.A = eye(2);
reset.b = 0;

guard = halfspace([1,0]',Text+18);

trans{1} = transition(guard,reset,1,'a','b'); % to on/on

guard = halfspace([0,-1]',-(Text+18));

trans{2} = transition(guard,reset,3,'a','b'); % to off/off

guard = mptPolytope([1,0;0,-1],[Text+18;Text+18]);

trans{3} = transition(guard,reset,4,'a','b'); % to on/off

loc{4} = location('on/off',4,inv,trans,linSys);

HA = hybridAutomaton(loc); % select location for hybrid automaton

% options
%set options---------------------------------------------------------------
options.x0 = [Text+16; Text+16]; %initial state for simulation
options.R0 = zonotope([options.x0, diag([1, 1])]); %initial set
options.startLoc = 1; %initial location
options.finalLoc = 0; %0: no final location
options.tStart = 0; %start time
options.tFinal = 24 * 3600; %final time
options.timeStepLoc{1} = 1; %time step size in location 1
options.timeStepLoc{2} = 1; %time step size in location 1
options.timeStepLoc{3} = 1; %time step size in location 1
options.timeStepLoc{4} = 1; %time step size in location 1
options.taylorTerms = 10; % reachability ?
options.polytopeType = 'mpt';

 options.zonotopeOrder = 20;
 options.polytopeOrder = 10;
 options.errorOrder=2;
 options.reductionTechnique = 'girard';
 options.isHybrid = 1;
 options.isHyperplaneMap = 0;
 options.enclosureEnables = [5]; %choose enclosure method(s)
 options.originContained = 0;
%--------------------------------------------------------------------------
%set input:
options.uLoc{1} = [1; 1]; %input for simulation
options.uLocTrans{1} = options.uLoc{1}; %center of input set
options.Uloc{1} = zonotope(zeros(2,1)); %input deviation from center

options.uLoc{2} = [1; 1]; %input for simulation
options.uLocTrans{2} = options.uLoc{2}; %center of input set
options.Uloc{2} = zonotope(zeros(2,1)); %input deviation from center

options.uLoc{3} = [1; 1]; %input for simulation
options.uLocTrans{3} = options.uLoc{3}; %center of input set
options.Uloc{3} = zonotope(zeros(2,1)); %input deviation from center

options.uLoc{4} = [1; 1]; %input for simulation
options.uLocTrans{4} = options.uLoc{4}; %center of input set
options.Uloc{4} = zonotope(zeros(2,1)); %input deviation from center

display('Begin simulation');

%simulate hybrid automaton
HA = simulate(HA,options);

display('End simulation');

% display('Begin reachability analysis');
% %compute reachable set
% [HA] = reach(HA,options);
% 
% display('End reachability analysis');
% 
%  %choose projection and plot------------------------------------------------
%  options.projectedDimensions = [1 2];
%  options.plotType = 'b';
%  plot(HA,'reachableSet',options); %plot reachable set
%  plot(options.R0,options.projectedDimensions,'blackFrame'); %plot initial set
%  plot(HA,'simulation',options); %plot simulation
%  %--------------------------------------------------------------------------

end