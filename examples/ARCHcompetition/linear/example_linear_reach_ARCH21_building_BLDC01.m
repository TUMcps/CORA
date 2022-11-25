function res = example_linear_reach_ARCH21_building_BLDC01()
% example_linear_reach_ARCH21_building_BLDC01 - building benchmark from
%    the 2021 ARCH competition (constant inputs)
%
% Syntax:  
%    res = example_linear_reach_ARCH21_building_BLDC01()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 

% Author:       Niklas Kochdumper
% Written:      13-March-2019
% Last update:  18-October-2019
% Last revision:---

%------------- BEGIN CODE --------------



% Parameters --------------------------------------------------------------

% input set
U = interval(0.8,1);

% initial set
X0 = interval([0.0002*ones(10,1); zeros(14,1); -0.0001; zeros(23,1)],....
              [0.00025*ones(10,1); zeros(14,1); 0.0001; zeros(23,1)]);
          
% model constant input as additional state
R0 = zonotope(cartProd(X0,U));


% Reachability Settings ---------------------------------------------------

options.taylorTerms = 4;
options.zonotopeOrder = 100;

options.linAlg = 'wrapping-free';

timeStep_1 = 0.002;
timeStep_2 = 0.01;


% System Dynamics ---------------------------------------------------------

% load system matrices 
load build

% define input as an additional state
dim = size(A,1);

% extended state matrix, input matrix all-zero
A_ = [[A,B];zeros(1,dim+1)];
B_ = zeros(dim+1,1);

% construct linear system objects
sys = linearSys('buildingSys',A_,B_);


% Reachability Analysis ---------------------------------------------------

timer = tic;

% first phase: time interval t \in [0,1] s
params.R0 = R0;
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

res = ~violation & boundReached1 & ~boundReached2;
tVer = toc(timer);

disp(['specifications verified: ',num2str(res)]);
disp(['computation time: ',num2str(tVer + tComp)]);



% Simulation --------------------------------------------------------------

% adapt parameter
params.R0 = R0;
params = rmfield(params,'tStart');

% settings for random simulation
simOpt.points = 10;
simOpt.fracVert = 0.0;
simOpt.fracInpVert = 0; % no inputs
simOpt.inpChanges = 0;  % no inputs

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


% Plot 2: time interval t \in [1,20] s
figure; hold on; box on;
plotOverTime(R,25,'FaceColor',[.6 .6 .6],'EdgeColor','none');
plotOverTime(simRes,25);
axis([0, 20, -8e-3, 6e-3])
xlabel('t');
ylabel('x_{25}');

close all;

%------------- END OF CODE --------------