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

inv = interval(temperatureInterval);

reset.A = eye(2);
reset.b = 0;

guard = halfspace([-1,0]',-(Text+22));

trans{1} = transition(guard,reset,4); % to off/on

guard = halfspace([0,-1]',-(Text+18));

trans{2} = transition(guard,reset,2); % to on/off

guard = mptPolytope([-1,0;0,-1],[-(Text+22);-(Text+18)]);

trans{3} = transition(guard,reset,3); % to off/off

loc{1} = location('on/on',inv,trans,linSys);

% Location 1B : Heater (on, off)
A = [-2*a,a;a,-2*a];
B = [a*Text + P1 / C,0;0,a*Text];

linSys = linearSys('lin',A,B);

inv = interval(temperatureInterval);

reset.A = eye(2);
reset.b = 0;

guard = halfspace([0,1]',Text+14);

trans{1} = transition(guard,reset,1); % to on/on

guard = halfspace([-1,0]',-(Text+22));

trans{2} = transition(guard,reset,3); % to off/off

guard = mptPolytope([-1,0;0,1],[-(Text+22);Text+14]);

trans{3} = transition(guard,reset,4); % to off/on

loc{2} = location('on/off',inv,trans,linSys);

% Location 2B : Heater (off,off)
A = [-2*a,a;a,-2*a];
B = [a*Text,0;0,a*Text];

linSys = linearSys('lin',A,B);

inv = interval(temperatureInterval);

reset.A = eye(2);
reset.b = 0;

guard = halfspace([1,0]',Text+18);

trans{1} = transition(guard,reset,2); % to on/off

guard = halfspace([0,1]',Text+14);

trans{2} = transition(guard,reset,4); % to off/on

guard = mptPolytope([1,0;0,1],[Text+18;Text+14]);

trans{3} = transition(guard,reset,1); % to on/on

loc{3} = location('on/off',inv,trans,linSys);

% Location 2A : Heater (off,on)
A = [-2*a,a;a,-2*a];
B = [a*Text,0;0,a*Text + P2 / C];

linSys = linearSys('lin',A,B);

inv = interval(temperatureInterval);

reset.A = eye(2);
reset.b = 0;

guard = halfspace([1,0]',Text+18);

trans{1} = transition(guard,reset,1); % to on/on

guard = halfspace([0,-1]',-(Text+18));

trans{2} = transition(guard,reset,3); % to off/off

guard = mptPolytope([1,0;0,-1],[Text+18;Text+18]);

trans{3} = transition(guard,reset,4); % to on/off

loc{4} = location('on/off',inv,trans,linSys);

HA = hybridAutomaton(loc); % select location for hybrid automaton

% options
%set options---------------------------------------------------------------
params.x0 = [Text+16; Text+16]; %initial state for simulation
params.R0 = zonotope([params.x0, diag([1, 1])]); %initial set
params.startLoc = 1; %initial location
params.finalLoc = 0; %0: no final location
params.tStart = 0; %start time
params.tFinal = 24 * 3600; %final time
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
options.enclose = [5]; %choose enclosure method(s)
options.originContained = 0;
%--------------------------------------------------------------------------
%set input:
params.uLoc{1} = [1; 1]; %input for simulation
params.uLocTrans{1} = params.uLoc{1}; %center of input set
params.Uloc{1} = zonotope(zeros(2,1)); %input deviation from center

params.uLoc{2} = [1; 1]; %input for simulation
params.uLocTrans{2} = params.uLoc{2}; %center of input set
params.Uloc{2} = zonotope(zeros(2,1)); %input deviation from center

params.uLoc{3} = [1; 1]; %input for simulation
params.uLocTrans{3} = params.uLoc{3}; %center of input set
params.Uloc{3} = zonotope(zeros(2,1)); %input deviation from center

params.uLoc{4} = [1; 1]; %input for simulation
params.uLocTrans{4} = params.uLoc{4}; %center of input set
params.Uloc{4} = zonotope(zeros(2,1)); %input deviation from center

display('Begin simulation');

%simulate hybrid automaton
HA = simulate(HA,params);

display('End simulation');

% display('Begin reachability analysis');
% %compute reachable set
% [HA] = reach(HA,params,options);
% 
% display('End reachability analysis');
% 
%  %choose projection and plot------------------------------------------------
%  options.projectedDimensions = [1 2];
%  options.plotType = 'b';
%  plot(HA,'reachableSet',options); %plot reachable set
%  plot(params.R0,options.projectedDimensions,'blackFrame'); %plot initial set
%  plot(HA,'simulation',options); %plot simulation
%  %--------------------------------------------------------------------------

end