function completed = example_linearParam_reach_01_rlc_const()
% example_linearParam_reach_01_rlc_const - example of linear parametric
%     reachability analysis with constant parameters, taken from [1]
%
% Syntax:
%    completed = example_linearParam_reach_01_rlc_const()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] M. Althoff, B. H. Krogh, and O. Stursberg. "Modeling, Design, and 
%        Simulation of Systems with Uncertainties", chapter Analyzing 
%        Reachability of Linear Dynamic Systems with Parametric 
%        Uncertainties, pages 69-94. Springer, 2011.

% Authors:       Matthias Althoff
% Written:       18-August-2016
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% System Dynamics ---------------------------------------------------------

% get matrix zonotopes of the model
[matZ_A,matZ_B] = RLCcircuit();
matI_A = intervalMatrix(matZ_A);
dim_x = dim(matZ_A,1);

% create linear parametric systems with constant parameters
sysMatZono  = linParamSys(matZ_A, eye(dim_x));
sysIntMat = linParamSys(matI_A, eye(dim_x));


% Parameters --------------------------------------------------------------

% compute initial set
u0 = intervalMatrix(0,0.2);             % range of voltages

intA = intervalMatrix(matZ_A);
invAmid = inv(center(intA.int));        % inverse of A

intB = intervalMatrix(matZ_B);
R0 = invAmid*intB*u0 + intervalMatrix(0,1e-3*ones(dim_x,1));

params.R0 = zonotope(interval(R0));     % convert initial set to zonotope

% uncertain inputs
u = intervalMatrix(1,0.01);
params.U = zonotope(interval(intB*u));

% final time
params.tFinal = 0.3;


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.001;
options.zonotopeOrder = 400;
options.taylorTerms = 8;
options.intermediateTerms = 2;
options.compTimePoint = false;


% Reachability Analysis ---------------------------------------------------

% compute reachable set using matrix zonotopes
tic
RmatZono = reach(sysMatZono, params, options);
tComp = toc;
disp(['computation time of reachable set using matrix zonotopes: ',num2str(tComp)]);

% compute reachable set using interval matrices
tic
RintMat = reach(sysIntMat, params, options);
tComp = toc;
disp(['computation time of reachable set using interval matrices: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 10;
simRes = simulateRandom(sysIntMat, params, simOpt);


% Visualization -----------------------------------------------------------

% PLOT 1: state space

figure; hold on; box on;
projDim = [20,40];
useCORAcolors("CORA:contDynamics", 2)
    
% plot reachable sets
plot(RintMat,projDim,'Order',10,'DisplayName','Interval matrix');
plot(RmatZono,projDim,'Order',10,'DisplayName','Matrix zonotope');

% plot initial set
plot(RintMat.R0,projDim,'DisplayName','Initial set');

% plot simulation results     
plot(simRes,projDim,'DisplayName','Simulations');

% label plot
xlabel(['x_{',num2str(projDim(1)),'}']);
ylabel(['x_{',num2str(projDim(2)),'}']);
legend();


% PLOT 2: reachable set over time

figure; hold on;
useCORAcolors("CORA:contDynamics", 2)

% plot time elapse
plotOverTime(RintMat,projDim(1),'DisplayName','Interval matrix');
plotOverTime(RmatZono,projDim(1),'DisplayName','Matrix zonotope');

% plot initial set
plotOverTime(RintMat.R0,projDim(1),'DisplayName','Initial set');

% plot simulation results
plotOverTime(simRes,projDim(1),'DisplayName','Simulations');

% label plot
xlabel('t');
ylabel(['x_{',num2str(projDim(1)),'}']);
legend();

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
