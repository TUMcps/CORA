function completed = example_linProbSys_reach_01_2dim()
% example_linProbSys_reach_01_2dim - example of probabilistic reachability 
%    analysis of a linear system with uncertain inputs, taken from 
%    [1, Sec. 4.2.8].
%
% Syntax:
%    completed = example_linProbSys_reach_01_2dim()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] M. Althoff, “Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, TUM 2010

% Authors:       Matthias Althoff
% Written:       16-July-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 2.5;
% probabilistic set of initial states
params.R0 = probZonotope([[3; 3],0.1*eye(2)],[0.7 0; 0 0.7],4);
params.U = zonotope([zeros(2,1),0.1*eye(2)]);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.01; 
options.taylorTerms = 4;
options.zonotopeOrder = 10; 
options.gamma = 4;


% System Dynamics ---------------------------------------------------------

%system matrix 
A = [-1 -4; 4 -1];
%input matrix 
B = eye(2);
%noise matrix 
C = 0.7*eye(2);

twoDimSys = linProbSys('twoDimSys',A,B,C);


% Reachability Analysis ---------------------------------------------------

tic;
R = reach(twoDimSys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 25;
simOpt.timeStep = options.timeStep;

simRes = simulateRandom(twoDimSys, params, simOpt);


% Visualization -----------------------------------------------------------

% plot different projections
dims = {[1 2]};

figure;

% plot reachable set
for k = 1:length(dims)
    
    subplot(1,2,1); hold on; box on
    projDims = dims{k};

    % plot reachable sets 
    plot(R,projDims,'FaceColor','next','m',2.5,'FaceColor','interp','DisplayName','Reachable set');
    
    % plot initial set
    plot(zonotope(params.R0,3),projDims,'k','FaceColor','w','ZPos',1,'DisplayName','Initial set');
    
    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
    legend('Location','east');

    % set view
    view(-35,30);
    
    % construct custom color map
    l=linspace(1,0,100)';
    colormap([l,l,l]);
    cmap = colormap;
    newCmap=[ones(1,3);cmap(2:end,:)];
    colormap(newCmap);
end


% plot simulation runs
for k = 1:length(dims)
    
    subplot(1,2,2); hold on; box on
    projDims = dims{k};
    
    % plot initial set
    plot(zonotope(params.R0,3),projDims,'k','FaceColor','w','ZPos',1,'DisplayName','Initial set');
    
    % plot simulation results
    plot(simRes,projDims,'k','DisplayName','Simulations');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
    legend('Location','east')

    % set view
    view(-35,30);
    xlim([-6,6]); ylim([-4,6]); zlim([0,1]);
end

% example completed
completed = true;


% ------------------------------ END OF CODE ------------------------------
