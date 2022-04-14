function completed = example_nonlinearParam_reach_02_tank_certainCase()
% example_nonlinearParam_reach_02_tank_certainCase - example of nonlinear
%    reachability analysis with uncertain parameters; this example is also 
%    a unit test function. The difference compared to the uncertain case is
%    that the parameter value is fixed. Unlike the example in the class
%    nonlinearSys, one can change the parameter values using the
%    options.paramInt.
%
% This example can be found in [1, Sec. 3.4.5] or in [2].
%
%
% Syntax:  
%    completed = example_nonlinearParam_reach_02_tank_certainCase()
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean 
%
% References:
%    [1] M. Althoff, “Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, TUM 2010
%    [2] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", CDC 2008 
% 
% Author:       Matthias Althoff
% Written:      19-August-2016
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---

%------------- BEGIN CODE --------------

% Parameters --------------------------------------------------------------

dim_x = 6;
params.tFinal = 400;
params.R0 = zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim_x)]); % initial set
params.U = zonotope([0,0.005]);                             % uncertain input


% Reachability Settings ---------------------------------------------------

options.timeStep = 1;
options.taylorTerms = 4;
options.zonotopeOrder = 10;

options.tensorOrder = 2;
options.alg = 'lin';


% System Dynamics ---------------------------------------------------------

% system without uncertain paramters 
tank = nonlinearSys(@tank6Eq);

% system with uncertain parameters
tankParam = nonlinParamSys(@tank6paramEq);


% Reachability Analysis ---------------------------------------------------
        
% compute reachable set of tank system without uncertain parameters
tic
RcontNoParam = reach(tank, params, options);
tComp = toc;
disp(['computation time of reachable set without uncertain parameters: ',num2str(tComp)]);

% compute reachable set of tank system with uncertain parameters
options.intermediateTerms = 4;
params.paramInt = 0.015;
tic
RcontParam = reach(tankParam, params, options);
tComp = toc;
disp(['computation time of reachable set with uncertain parameters: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 60;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;
params = rmfield(params,'paramInt');

simRes = simulateRandom(tank, params, simOpt);


% Visualization -----------------------------------------------------------

% plot different projections
dims = {[1 2],[3 4],[5 6]};

for k = 1:length(dims)

    figure;
    projDims = dims{k};
    
    subplot(1,2,1); hold on; box on

    % plot reachable set (with uncertain paramaeters)
    plot(RcontParam,projDims,'FaceColor',[.6 .6 .6],'EdgeColor','none');

    % plot initial set
    plot(params.R0,projDims,'w','Filled',true,'EdgeColor','k');
    
    subplot(1,2,2); hold on; box on
    
    % plot reachable set (without uncertain parameters)
    plot(RcontNoParam,projDims,'w','EdgeColor','k');

    % plot initial set
    plot(params.R0,projDims,'w','Filled',true,'EdgeColor','k');
    
    % loop over all subplots
    for sp=1:2
        subplot(1,2,sp); hold on;
        
        % plot simulation results      
        plot(simRes,projDims,'b');

        % label plot
        xlabel(['x_{',num2str(projDims(1)),'}']);
        ylabel(['x_{',num2str(projDims(2)),'}']);
    end
end

% example completed
completed = 1;

%------------- END OF CODE --------------
