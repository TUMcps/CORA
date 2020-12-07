function completed = example_linearDT_reach_01_5dim()
% example_linearDT_reach_01_5dim - example of discrete-time linear 
%                                  reachability analysis with uncertain 
%                                  inputs
%
% Syntax:  
%    example_linearDT_reach_01_5dim
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
% 
% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      20-March-2020
% Last update:  23-April-2020 (restructure params/options)
%               07-December-2020 (add output matrix)
% Last revision:---


%------------- BEGIN CODE --------------


% System Dynamics ---------------------------------------------------------

% system matrix
A = [0.95 -0.15 0 0 0;...
    0.15 0.95 0 0 0; ...
    0 0 0.9 0.05 0; ...
    0 0 -0.05 0.9 0; ...
    0 0 0 0 0.92];

% input matrix
B = 1;

% constant input
c = zeros(5,1);

% output matrix
C = [2 0 0 0 0; 0 1 0 0 0];

% sampling time
dt = 0.04;

fiveDimSys = linearSysDT('fiveDimSys',A,B,c,C,dt);


% Parameter ---------------------------------------------------------------

params.tFinal = 5;
params.R0 = zonotope([ones(length(A),1),0.1*eye(length(A))]);
params.U = zonotope([zeros(5,1),0.02*diag([0.1, 0.3, 0.1, 0.3, 0.3])]);


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 200;


% Reachability Analysis ---------------------------------------------------

tic
R = reach(fiveDimSys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 25;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 10;

simRes = simulateRandom(fiveDimSys, params, simOpt);

% apply output equation y = Cx
for i=1:simOpt.points
    outputtraj{i,1} = (C*simRes.x{i}')';
end
simRes = simResult(outputtraj,simRes.t);

% Visualization -----------------------------------------------------------

% plot different projections
dims = {[1 2]};

for k = 1:length(dims)
    
    figure; hold on; box on
    projDims = dims{k};

    % plot reachable set
    plot(R,projDims,'FaceColor',[.8 .8 .8],'EdgeColor','b');
    
    % plot initial set
    plot(C*params.R0,projDims,'w-','lineWidth',2);
    
    % plot simulation results
    plot(simRes,projDims,'.k');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
completed = 1;

%------------- END OF CODE --------------