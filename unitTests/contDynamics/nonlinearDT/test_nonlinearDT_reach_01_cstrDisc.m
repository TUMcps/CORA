function res = test_nonlinearDT_reach_01_cstrDisc(~)
% test_nonlinearDT_reach_cstrDisc - unit-test for nonlinear discrete time
% reachability analysis.
%
% This example can be found in Sec. 6 of
% J.M. Bravo, Robust MPC of constrained discrete-time nonlinear systems 
% based on approximated reachable sets, 2006
%
% Syntax:  
%    test_nonlinearDT_reach_cstr
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

dim = 2;

% set options -------------------------------------------------------------
options.tStart=0;           % start time
options.tFinal=0.03;        % final time
options.x0=[-0.15;-45];     % initial state for simulation
options.R0=zonotope([options.x0,diag([0.005;3])]); % initial state for reachability analysis

options.taylorTerms=20;                 % number of taylor terms for reachable sets
options.zonotopeOrder=100;              % zonotope order in reachability analysis
options.intermediateOrder=50;           % intermediate zonotope order used in reachability analysis
options.reachabilitySteps=3;           % number of reachability steps in each intermediate time step
options.timeStep=options.tFinal...
        / options.reachabilitySteps;    % length of each time step in reachability analysis

% additional parameters for reachability analysis
options.advancedLinErrorComp = 0;
options.tensorOrder = 3;
options.errorOrder = 5;
options.reductionInterval = inf;
options.loc=1;
options.reductionTechnique='girard';
options.maxError = 150e5*ones(dim,1);
options.oldError = zeros(dim,1);

options.originContained = 0;


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
runs = 50;

simRes = simulate_random(sysDisc, options, runs, fractionVertices, fractionInputVertices);


% verify reachable set with simulation ------------------------------------
res = 1;

for i = 1:runs
   for j = 1:length(R)
      if ~containsPoint(R{j},simRes.x{i}(:,j))
          res = 0;
      end
   end
end

%------------- END OF CODE --------------