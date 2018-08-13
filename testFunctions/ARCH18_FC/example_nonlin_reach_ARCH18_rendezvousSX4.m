function example_nonlin_reach_ARCH18_rendezvousSX4()
% example_nonlin_reach_ARCH18_rendezvousSX4 - example of nonlinear 
% reachability analysis from the ARCH18 friendly competition 
% (instance of rendevouz example)
%
% Syntax:  
%    example_nonlin_reach_ARCH18_rendezvousSX4
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
% Written:      26-April-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% x_1 = v_x
% x_2 = v_y
% x_3 = s_x 
% x_4 = s_y

% debugging mode
options.debug = 0;

% initial set
R0 = zonotope([[-65; -25; 2; 0.8; 0],diag([5,5,0.05,0.05,0])]);
options.x0=center(R0); 
options.R0=R0;

% system input
options.U = zonotope(0);
options.uTrans = 0;
options.u = 0;

% parameters for reachability analysis
options.tStart=0;                   % start time
options.tFinal=100;                 % duration of reachability analysis
options.taylorTerms=20;             % number of taylor terms for reachable sets
options.zonotopeOrder=100;          % zonotope order in reachability analysis
options.intermediateOrder=50;       % intermediate zonotope order used in reachability analysis
options.reachabilitySteps=100;      % number of reachability steps 
options.timeStep=options.tFinal/options.reachabilitySteps; % length of each time step in reachability analysis

options.advancedLinErrorComp = 0;
options.tensorOrder = 2;
options.errorOrder = 5;
options.reductionInterval = inf;
options.loc=1;
options.reductionTechnique='girard';
options.maxError = 150e5*ones(5,1);
options.oldError = zeros(5,1);

% load system
dynamics = nonlinearSys(5,1,@rendezvous_nonlin_4d,options);

% reachability analysis
[Rcont,R] = reach(dynamics,options);

% simulation
[~,~,xSim] = simulate(dynamics,options,options.tStart,options.tFinal,options.x0,options);

% plot the result
hold on
for i = 1:length(Rcont)
   plotFilled(Rcont{i}{1},[1,2],'g'); 
end
plot(xSim(:,1),xSim(:,2),'k');
fill([-0.1,0.1,0.1,-0.1],[-0.1,-0.1,0.1,0.1],'r','EdgeColor','none');
xlabel('x [m]');
ylabel('y [m]');
grid on

figure
hold on
for i = 1:length(R)
   plotFilled(Rcont{i}{1},[3,4],'g'); 
end
plot(xSim(:,3),xSim(:,4),'k');
xlabel('v_x [m/min]');
ylabel('v_y [m/min]');
grid on

%------------- END OF CODE --------------
