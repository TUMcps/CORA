function res = example_linear_reach_ARCH21_building_BLDF01()
% example_linear_reach_ARCH21_building_BLDF01 - building benchmark from the
%    ARCH 2021 competition (time-varying inputs)
%
% Syntax:  
%    res = example_linear_reach_ARCH21_building_BLDF01()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%

% Author:       Niklas Kochdumper
% Written:      13-March-2019
% Last update:  18-October-2019
% Last revision:---

%------------- BEGIN CODE --------------


% Parameters --------------------------------------------------------------

% initial set 
R0 = interval([0.0002*ones(10,1); zeros(14,1); -0.0001; zeros(23,1)],...
              [0.00025*ones(10,1); zeros(14,1); 0.0001; zeros(23,1)]);
      
% uncertain inputs
params.U = zonotope(interval(0.8,1));


% Reachability Settings ---------------------------------------------------

options.taylorTerms = 4;
options.zonotopeOrder = 100;

options.linAlg = 'wrapping-free';

timeStep_1 = 0.002;
timeStep_2 = 0.01;


% System Dynamics ---------------------------------------------------------

% load system matrices
load build

% construct linear system object
sys = linearSys('buildingSys',A,B);


% Reachability Analysis ---------------------------------------------------

timer = tic;

% first phase: time interval t \in [0,1] s
params.R0 = zonotope(R0);
params.tFinal = 1;

options.timeStep = timeStep_1;

R1 = reach(sys, params, options);

% second phase: time interval t \in [1,20] s
params.R0 = R1.timePoint.set{end};
params.tStart = 1;
params.tFinal = 20;

options.timeStep = timeStep_2;

R2 = reach(sys, params, options);

tComp = toc(timer);

R = add(R1,R2);


% Verification ------------------------------------------------------------

timer = tic;
violation = 0;
boundReached1 = 0;
boundReached2 = 0;

for i=1:length(R1.timeInterval.set)
    x_25 = interval(project(R1.timeInterval.set{i},25));
    if supremum(x_25) > 5.1e-3
        violation = 1;
    end
    if supremum(x_25) > 4e-3
        boundReached1 = 1;
    end
end
for i=1:length(R2.timeInterval.set)
    x_25 = interval(project(R2.timeInterval.set{i},25));
    if supremum(x_25) > 5.1e-3
        violation = 1;
    end
end
x_25 = interval(project(R2.timeInterval.set{end},25));
if infimum(x_25) < -0.78e-3
   boundReached2 = 1; 
end

res = ~violation & boundReached1 & boundReached2;
tVer = toc(timer);

disp(['specifications verified: ',num2str(res)]);
disp(['computation time: ',num2str(tVer + tComp)]);



% Simulation --------------------------------------------------------------

% adapt parameter
params.R0 = zonotope(R0);
params = rmfield(params,'tStart');

% settings for random simulation
simOpt.points = 10;
simOpt.fracVert = 0.0;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 20;

% simulate the system
simRes = simulateRandom(sys, params, simOpt);



% Visualization -----------------------------------------------------------

% Plot 1: time interval t \in [0,1] s
figure; hold on; box on;
plotOverTime(R1,25,'FaceColor',[.6 .6 .6],'EdgeColor','none');
plotOverTime(simRes,25);
axis([0, 1, -8e-3, 6e-3])
xlabel('t');
ylabel('x_{25}');
box on


% Plot 2: time interval t \in [0,20] s
figure; hold on; box on;
plotOverTime(R,25,'FaceColor',[.6 .6 .6],'EdgeColor','none');
plotOverTime(simRes,25);
axis([0, 20, -8e-3, 6e-3])
xlabel('t');
ylabel('x_{25}');
box on

%------------- END OF CODE --------------