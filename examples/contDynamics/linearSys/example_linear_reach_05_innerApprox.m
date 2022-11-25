function res = example_linear_reach_05_innerApprox()
% example_linear_reach_05_innerApprox - example that demonstrates how to
%                                       compute an inner-approximation of 
%                                       the reachable set
%
% Syntax:  
%    res = example_linear_reach_05_innerApprox()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Author:       Niklas Kochdumper
% Written:      27-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% Parameter ---------------------------------------------------------------

params.tFinal = 1;
params.R0 = zonotope([ones(5,1),0.1*diag(ones(5,1))]);
params.U = zonotope(interval([0.9; -0.25; -0.1; 0.25; -0.75], ...
                             [1.1; 0.25; 0.1; 0.75; -0.25]));


% Reachability Settings ---------------------------------------------------

% settings for inner-approximation
options1.timeStep = 0.02;
options1.zonotopeOrder = 20;

% settings for outer-approximation
options2.timeStep = 0.02; 
options2.taylorTerms = 4;
options2.zonotopeOrder = 20; 


% System Dynamics ---------------------------------------------------------

A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;

fiveDimSys = linearSys('fiveDimSys',A,B);


% Reachability Analysis ---------------------------------------------------

% compute inner-approximation
tic
Rin = reachInner(fiveDimSys, params, options1);
tComp = toc;

disp(['computation time (inner-approximation): ',num2str(tComp)]);

% compute outer-approximation
tic
Rout = reach(fiveDimSys, params, options2);
tComp = toc;

disp(['computation time (outer-approximation): ',num2str(tComp)]);



% Visualization -----------------------------------------------------------

% plot different projections
dims = {[1 2],[3 4]};

% Plot 1: reachable set for whole time horizon
for k = 1:length(dims)

    figure; hold on; box on
    projDims = dims{k};

    % plot outer-approximation 
    for i = 1:length(Rout.timePoint.set)
        h1 = plot(Rout.timePoint.set{i},projDims,'b');
    end
    
    % plot inner-approximation
    h2 = plot(Rin,projDims,'r','Filled',false);

    % plot initial set
    plot(params.R0,projDims,'w','Filled',true);

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
    
    % legend
    legend([h1,h2],'outer-approx.','inner-approx.');
end

% Plot 2: final reachable set
for k = 1:length(dims)

    figure; hold on; box on
    projDims = dims{k};

    % plot outer-approximation 
    h1 = plot(Rout.timePoint.set{end},projDims,'b');
    
    % plot inner-approximation
    h2 = plot(Rin.timePoint.set{end},projDims,'r','Filled',false);

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
    
    % legend
    legend([h1,h2],'outer-approx.','inner-approx.');
end

% example completed
res = 1;

%------------- END OF CODE --------------