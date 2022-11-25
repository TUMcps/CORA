function completed = example_linearDT_reach_02_5dim()
% example_linearDT_reach_02_5dim - example for linear reachability 
%                                  analysis using discrete time
%
% Syntax:  
%    completed = example_linearDT_reach_02_5dim()
%
% Inputs:
%    no
%
% Outputs:
%    completed - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      31-July-2018
% Last update:  23-April-2020 (moved here from contDynamics/linearSys)
% Last revision:---


%------------- BEGIN CODE --------------

% System Dynamics ---------------------------------------------------------

% system matrices
A = [0.8, -0.4, 0,    0,    0; ...
     0.4,  0.8, 0,    0,    0; ...
     0,    0,   0.7,  0.07, 0; ...
     0,    0,  -0.07, 0.7,  0; ...
     0,    0,   0,    0,    0.8];
dim_x = length(A);
 
B = eye(dim_x);

% sampling time
dt = 0.04;

% discrete-time linear system object
sys = linearSysDT('fiveDimSys',A,B,dt);


% Parameter ---------------------------------------------------------------

params.tFinal = 5;
params.R0 = zonotope(ones(dim_x,1),0.1*eye(dim_x));
params.U = zonotope(zeros(dim_x,1), 0.1*diag([0.2, 0.5, 0.2, 0.5, 0.5]));


% Reachability Settings  --------------------------------------------------

options.zonotopeOrder = 10;


% Reachability Analysis ---------------------------------------------------

tic
R = reach(sys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 25;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 10;

simRes = simulateRandom(sys, params, simOpt);


% Visualization -----------------------------------------------------------

% plot different projections
dims = {[1 2],[3 4]};

for k = 1:length(dims)
    
    figure; hold on; box on
    projDims = dims{k};

    % plot reachable set
    plot(R,projDims,'FaceColor',[.8 .8 .8],'EdgeColor','b');
    
    % plot initial set
    plot(params.R0,projDims,'w-','lineWidth',2);
    
    % plot simulation results
    plot(simRes,projDims,'.k');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
completed = 1;

%------------- END OF CODE --------------
