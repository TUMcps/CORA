function completed = example_nonlinearDT_reach_cstrDisc
% example_nonlinearDT_reach_cstrDisc - example of nonlinear discrete time
% reachability analysis.
%
% This example can be found in Sec. 6 of
% J.M. Bravo, Robust MPC of constrained discrete-time nonlinear systems 
% based on approximated reachable sets, 2006
%
% Syntax:  
%    example_nonlinearDT_reach_cstr
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Niklas Kochdumper
% Written:      30-January-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


% set options -------------------------------------------------------------
options.tStart=0;           % start time
options.tFinal=0.15;        % final time
options.x0=[-0.15;-45];     % center point of the initial set
options.R0=zonotope([options.x0,diag([0.005;3])]); % initial set 

options.zonotopeOrder=100;              % maximum zonotope order
options.reachabilitySteps=10;           % number of reachability steps 
options.timeStep=options.tFinal...
        / options.reachabilitySteps; 

% additional parameters for reachability analysis
options.tensorOrder = 3;
options.errorOrder = 5;
options.reductionTechnique='girard';


% obtain uncertain inputs -------------------------------------------------
options.uTrans = [0;0];
options.U = zonotope([zeros(2,1),diag([0.1;2])]); %input for reachability analysis


% specify discrete dynamics ----------------------------------------------- 
sysDisc = nonlinearSysDT(2,2,@cstrDiscr,options);


% compute reachable sets --------------------------------------------------
R = reach(sysDisc,options);


% simulate the system -----------------------------------------------------
fractionVertices = 0.5;
fractionInputVertices = 0.5;
runs = 100;

simRes = simulate_random(sysDisc, options, runs, fractionVertices, fractionInputVertices);


% plot reachable sets -----------------------------------------------------
hold on
for i=1:length(R)
    Zproj = reduce(R{i},'girard',3);
    plotFilled(Zproj,[1 2],[.8 .8 .8],'EdgeColor','none');
end

% plot simulated trajectories ---------------------------------------------
for i=1:length(simRes.x)
    plot(simRes.x{i}(1,:),simRes.x{i}(2,:),'.k');
end

% add labels to plot ------------------------------------------------------
xlabel('T-T_0');
ylabel('C-C_0');
box on

completed = 1;

%------------- END OF CODE --------------