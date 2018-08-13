function PHA = testParallelHybridAutomaton()

% Test parallelHybridAutomaton

% Constant
a = 7.9*10^-5;
Text = 273.15; % 0K
% Text = 0;
P = 210*10^3;
C = 7.6*10^7;

% Room 1
% Location 1 : Heating on
A = [-2*a];
B = [a,a*Text + P / C];

linSys = linearSys('lin',A,B);

inv = intervalhull([Text+0,Text+100]);

%toff = 22C
guard = halfspace(-1,-(Text+22));

reset.A = 1;
reset.b = 0;

trans = {transition(guard,reset,2,'a','b')}; %transition

loc{1} = location('on',1,inv,trans,linSys);

% Location 2 : Heating off
A = [-2*a];
B = [a,a*Text];

linSys = linearSys('lin',A,B);

inv = intervalhull([Text+0,Text+100]);

%ton = 18C
guard = halfspace(1,Text+18);

reset.A = 1;
reset.b = 0;

trans = {transition(guard,reset,1,'a','b')}; %transition

loc{2} = location('off',2,inv,trans,linSys);

bd = bind(1,2,1);

comp{1} = component('Room 1',[], loc, {bd});

% Room 2
P = 160*10^3;
% Location 1 : Heating on
A = [-2*a];
B = [a,a*Text + P / C];

linSys = linearSys('lin',A,B);

inv = intervalhull([Text+0,Text+100]);
% Toff 18C
guard = halfspace(-1,-(Text+18));

reset.A = 1;
reset.b = 0;

trans = {transition(guard,reset,2,'a','b')}; %transition

loc{1} = location('on',1,inv,trans,linSys);

% Location 2 : Heating off
A = [-2*a];
B = [a,a*Text];

linSys = linearSys('lin',A,B);

inv = intervalhull([Text+0,Text+100]);

%Ton = 14C
guard = halfspace(1,Text+14);

reset.A = 1;
reset.b = 0;

trans = {transition(guard,reset,1,'a','b')}; %transition

loc{2} = location('off',2,inv,trans,linSys);

bd = bind(1,1,1);

comp{2} = component('Room 2',[], loc, {bd});

% Run
PHA = parallelHybridAutomaton(comp);

% options
%set options---------------------------------------------------------------
options.x0 = [Text+16; Text+16]; %initial state for simulation
options.R0 = zonotope([options.x0, diag([1, 1])]); %initial set
options.startLoc = [1 1]; %initial location
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
PHA = simulate(PHA,options);

display('End simulation');

% display('Begin reachability analysis');
% 
% %compute reachable set
% [PHA] = reach(PHA,options);

end