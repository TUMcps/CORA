function completed = example_nonlinear_reach_04_laubLoomis_polyZonotope()
% example_nonlinear_reach_04_laubLoomis_polyZonotope - example of nonlinear
%                                                      reachability analysis
% 
%    Example from [1] demonstrating reachability analysis for nonlinear
%    systems using polynomial zonotopes for the Laub-Loomis system.
%
% Syntax:  
%    completed = example_nonlinear_reach_04_laubLoomis_polyZonotope()
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean 
%
% References: 
%   [1] N. Kochdumper et al. "Sparse Polynomial Zonotopes: A Novel Set 
%       Representation for Reachability Analysis"

% Author:       Matthias Althoff
% Written:      18-August-2016
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---

%------------- BEGIN CODE --------------

% Parameters --------------------------------------------------------------

x0 = [1.2; 1.05; 1.5; 2.4; 1; 0.1; 0.45];
R0 = zonotope([x0,0.15*eye(7)]);
params.R0 = polyZonotope(R0);                          % initial set
params.tFinal = 20;                                    % final time


% Reachability Settings ---------------------------------------------------

% settings
options.timeStep = 0.01;
options.taylorTerms = 20;
options.zonotopeOrder = 100;
options.intermediateOrder = 50;
options.errorOrder = 15;

% reachability algorithm
options.alg = 'poly';
options.tensorOrder = 3;

% settings for polynomial zonotopes
polyZono.maxDepGenOrder = 30;
polyZono.maxPolyZonoRatio = 0.001;
polyZono.restructureTechnique = 'reduceFullGirard';

options.polyZono = polyZono;


% System Dynamics ---------------------------------------------------------

sys = nonlinearSys(@laubLoomis);



% Reachability Analysis ---------------------------------------------------

tic
R = reach(sys, params, options);
tComp = toc;

disp(['computation time of reachable set: ',num2str(tComp)]);



% Simulation --------------------------------------------------------------

% number of initial points
simOpt.points = 60;
% random simulation
simRes = simulateRandom(sys,params,simOpt); 



% Visualization -----------------------------------------------------------

% PLOT 1: reachable set over time

figure; hold on; box on

% plot reachable set
plotOverTime(R,4,'FaceColor',[0 0.8 0]);

% plot the unsafe set
I = interval([0;5],[params.tFinal;10]);
plot(I,[1,2],'FaceColor',[1 0.5 0],'FaceAlpha',0.6);

% plot simulation
plotOverTime(simRes,4);

% formatting
xlabel('t [s]');
ylabel('x_4');
xlim([0,params.tFinal]);
ylim([1,6]);


% PLOTs 2-3: projections of the reachable set

dims = {[1,2],[3,5],[6,7]};

for i = 1:length(dims)
    
    figure; hold on; box on;
    projDim = dims{i};
    
    % plot reachable set
    plot(R,projDim,'FaceColor',[0 0.8 0]);
    
    % plot initial set
    plot(R0,projDim,'k','FaceColor','w');
    
    % plot simulation results
    plot(simRes,projDim);
    
    % formatting
    xlabel(['x_',num2str(projDim(1))]);
    ylabel(['x_',num2str(projDim(2))]);
end

completed = true;

%------------- END OF CODE --------------