function completed = example_linearParam_reach_01_rlc_const()
% example_linearParam_reach_01_rlc_const - example of linear parametric
%                                          reachability analysis with 
%                                          constant parameter
%
% This example is taken from [1].
%
% Syntax:  
%    completed = example_linearParam_reach_01_rlc_const()
%
% Inputs:
%    no
%
% Outputs:
%    completed - boolean 
%
% References:
%    [1] M. Althoff, B. H. Krogh, and O. Stursberg. "Modeling, Design, and 
%        Simulation of Systems with Uncertainties", chapter Analyzing 
%        Reachability of Linear Dynamic Systems with Parametric 
%        Uncertainties, pages 69-94. Springer, 2011.
% 
% Author:       Matthias Althoff
% Written:      18-August-2016
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---

%------------- BEGIN CODE --------------


% System Dynamics ---------------------------------------------------------

% get matrix zonotopes of the model
[matZ_A,matZ_B] = RLCcircuit();
matI_A = intervalMatrix(matZ_A);
dim_x = matZ_A.dim;

% create linear parametric systems with constant parameters
sysMatZono  = linParamSys(matZ_A, eye(dim_x));
sysIntMat = linParamSys(matI_A, eye(dim_x));


% Parameter ---------------------------------------------------------------

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

simOpt.points = 60;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;

simRes = simulateRandom(sysIntMat, params, simOpt);


% Visualization -----------------------------------------------------------

% PLOT 1: state space

figure; hold on; box on;
projDim = [20,40];
    
% plot reachable sets
hanIntMat = plot(RintMat,projDim,'FaceColor',[.6 .6 .6],'EdgeColor','none','Order',10);
hanMatZono = plot(RmatZono,projDim,'FaceColor',[.8 .8 .8],'EdgeColor','none','Order',10);

% plot initial set
plot(params.R0,projDim,'w','Filled',true,'EdgeColor','k');

% plot simulation results     
plot(simRes,projDim);

% label plot
xlabel(['x_{',num2str(projDim(1)),'}']);
ylabel(['x_{',num2str(projDim(2)),'}']);
legend([hanIntMat,hanMatZono],'Interval matrix','Matrix zonotope');


% PLOT 2: reachable set over time

figure; hold on;

% plot time elapse
hanIntMat = plotOverTime(RintMat,0.5*dim_x,'FaceColor',[.6 .6 .6],'EdgeColor','none');
hanMatZono = plotOverTime(RmatZono,0.5*dim_x,'FaceColor',[.8 .8 .8],'EdgeColor','none');

% plot simulation results
plotOverTime(simRes,0.5*dim_x);

% label plot
xlabel('t');
ylabel(['x_{',num2str(0.5*dim_x),'}']);
legend([hanIntMat,hanMatZono],'Interval matrix','Matrix zonotope');

% example completed
completed = 1;

%------------- END OF CODE --------------