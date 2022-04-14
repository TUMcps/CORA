function res = example_linear_reach_01_5dim()
% example_linear_reach_01_5dim - example of linear reachability 
%                                analysis with uncertain inputs
%
% This example can be found in [1, Sec. 3.2.3].
%
% Syntax:  
%    completed = example_linear_reach_01_5dim()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% References:
%    [1] M. Althoff, “Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, TUM 2010
% 
% Author:       Matthias Althoff
% Written:      17-August-2016
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---

%------------- BEGIN CODE --------------

% Parameter ---------------------------------------------------------------

params.tFinal = 5;
params.R0 = zonotope([ones(5,1),0.1*diag(ones(5,1))]);
params.U = zonotope(interval([0.9; -0.25; -0.1; 0.25; -0.75], ...
                             [1.1; 0.25; 0.1; 0.75; -0.25]));


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.02; 
options.taylorTerms = 4;
options.zonotopeOrder = 20; 


% System Dynamics ---------------------------------------------------------

A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;

fiveDimSys = linearSys('fiveDimSys',A,B);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(fiveDimSys, params, options);
tComp = toc;

disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 25;
simOpt.inpChanges = 10;
simOpt.p_conf = 0.8;

simRes = simulateRandom(fiveDimSys, params, simOpt, 'gaussian');


% Visualization -----------------------------------------------------------

% plot different projections
dims = {[1 2],[3 4]};

for k = 1:length(dims)

    figure; hold on; box on
    projDims = dims{k};

    % plot reachable sets 
    plot(R,projDims,'FaceColor',[.8 .8 .8],'EdgeColor','b');

    % plot initial set
    plot(params.R0,projDims,'w-','lineWidth',2);

    % plot simulation results
    plot(simRes,projDims,'y');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
res = 1;

%------------- END OF CODE --------------