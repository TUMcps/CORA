function completed = example_linProbSys_reach_02_5dim()
% example_linProbSys_reach_02_5dim - example of probabilistic reachability 
%    analysis of a linear system with uncertain inputs (five dimensions),
%    taken from [1, Sec. 4.2.8].
%
% Syntax:
%    completed = example_linProbSys_reach_02_5dim()
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

params.tFinal = 2;
% probabilistic set of initial states
params.R0 = probZonotope([[2;2;2;2;2],0.1*eye(5)],0.3*eye(5),3);
params.U = zonotope([[1; 0; 0; 0.5; -0.5],diag([0.2, 0.5, 0.2, 0.5, 0.5])]);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.04; 
options.taylorTerms = 4;
options.zonotopeOrder = 10; 
options.gamma = 3;


% System Dynamics ---------------------------------------------------------

% system matrix 
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
% input matrix 
B = eye(5);
% noise matrix 
C = 0.4*eye(5);

fiveDimSys = linProbSys('fiveDimSys',A,B,C);


% Reachability Analysis ---------------------------------------------------

tic;
R = reach(fiveDimSys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 25;
simOpt.timeStep = options.timeStep;

simRes = simulateRandom(fiveDimSys, params, simOpt);


% Visualization -----------------------------------------------------------

% plot different projections
dims = {[2 3], [4 5]};

figure

% plot reachable set
for k = 1:length(dims)
    
    subplot(2,2,k); hold on; box on
    projDims = dims{k};

    % plot reachable sets 
    plot(R,projDims,'Color','next','m',2.5,'FaceColor','interp','DisplayName','Reachable set');
    
    % plot initial set
    plot(zonotope(params.R0,3),projDims,'k','FaceColor','w','ZPos',5,'DisplayName','Initial set');
    
    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
    legend('Location','east')

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
    
    subplot(2,2,2+k); hold on; box on
    projDims = dims{k};
    
    % plot initial set
    plot(zonotope(params.R0,3),projDims,'k','FaceColor','w','ZPos',5,'DisplayName','Initial set');
    
    % plot simulation results
    plot(simRes,projDims,'k','DisplayName','Simulations');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
    legend('Location','east')

    % set view
    view(-35,30);
end

% example completed
completed = true;


% ------------------------------ END OF CODE ------------------------------
