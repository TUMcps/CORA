function completed = example_nonlinearDT_reach_cstrDisc
% example_nonlinearDT_reach_cstrDisc - example of nonlinear discrete time
% reachability analysis.
%
% This example can be found in [1, Sec. 6].
% 
%
% Syntax:  
%    example_nonlinearDT_reach_cstr
%
% Inputs:
%    no
%
% Outputs:
%    completed - boolean 
%
% Example: 
%
% References:
%    [1] J.M. Bravo, Robust MPC of constrained discrete-time
%        nonlinear systems based on approximated reachable sets, 2006
% 
% Author:       Niklas Kochdumper, Matthias Althoff
% Written:      30-January-2018
% Last update:  20-March-2020 (MA, simulateRandomDT from inherited class)
%               23-April-2020 (restructure params/options)
% Last revision:---


%------------- BEGIN CODE --------------

% Parameter  --------------------------------------------------------------

params.tFinal = 0.15;
params.U = zonotope([zeros(2,1),diag([0.1;2])]);
% zonotope
params.R0 = zonotope([[-0.15;-45],diag([0.005;3])]);
% zonotope bundle
% params.R0 = zonoBundle({ zonotope([[-0.15;-45],diag([0.005;5])]); ...
%     zonotope([[-0.15;-45],diag([0.005;3])]) });
% constrained zonotope
params.R0 = conZonotope([[-0.15;-45],[0.005 0 0; 0 1.5 1.5]], [1 1 -1], 1);


% Reachability Settings  --------------------------------------------------

options.zonotopeOrder = 10;
options.tensorOrder = 2;
options.errorOrder = 5;

% System Dynamics  --------------------------------------------------------

% sampling time
dt = 0.015;

fun = @(x,u) cstrDiscr(x,u,dt);

sysDisc = nonlinearSysDT('stirredTankReactor',fun,0.015);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(sysDisc,params,options);
tComp = toc;
disp("Computation time: " + tComp);


% Simulation --------------------------------------------------------------

simOpt.points = 10;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 3;

simRes = simulateRandom(sysDisc, params, simOpt);


% Visualization -----------------------------------------------------------

figure; hold on; box on;

% plot initial set
plot(params.R0,[1,2],'w','Filled',true);

% plot reachable set
if isa(params.R0,'conZonotope')
    plot(R,[1 2],'b','EdgeColor','none','Splits',5); % for conZonotope
else
    plot(R,[1 2],'b','EdgeColor','none');
end

% plot simulation
plot(simRes,[1,2],'.k');

% formatting
xlabel('T-T_0');
ylabel('C-C_0');


% example completed
completed = 1;

%------------- END OF CODE --------------