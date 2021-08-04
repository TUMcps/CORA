function completed = example_nonlinearParam_reach_01_tank()
% example_nonlinearParam_reach_01_tank - example of nonlinear reachability 
%                                        analysis with uncertain parameters
%
% This file implements the example from [1], which can also be found in 
% Sec. 3.4.5 of [2].
%
% Syntax:  
%    example_nonlinearParam_reach_01_tank
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% References:
%   [1] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%       uncertain parameters using conservative linearization"  
%   [2] M. Althoff “Reachability analysis and its application to the safety 
%       assessment of autonomous cars"

% Author:       Matthias Althoff
% Written:      19-August-2016
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---


%------------- BEGIN CODE --------------

% Parameter ---------------------------------------------------------------

params.tFinal = 400;                                     % final time
params.R0 = zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(6)]);  % initial set
params.U = zonotope([0,0.005]);                          % uncertain input

params_ = params;
params_.paramInt = interval(0.0148,0.015);          % uncertain paramters


% Reachability Settings ---------------------------------------------------

options.timeStep=0.5;
options.taylorTerms=4;
options.zonotopeOrder=10;
options.tensorOrder = 2;
options.alg = 'lin';


% System Dynamics ---------------------------------------------------------

% tank system with certain pararmters
tank = nonlinearSys(@tank6Eq);

% tank system with uncertain parameters
tankParam = nonlinParamSys(@tank6paramEq);


% Reachability Analysis ---------------------------------------------------        

% compute reachable set of tank system without uncertain parameters
tic
RcontNoParam = reach(tank, params, options);
tComp = toc;
disp(['computation time of reachable set without uncertain parameters: ',num2str(tComp)]);

% compute reachable set of tank system with uncertain parameters
options.intermediateTerms = 4;
tic
RcontParam = reach(tankParam, params_, options);
tComp = toc;
disp(['computation time of reachable set with uncertain parameters: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

% settings for random simulation 
simOpt.points = 60;        % number of initial points
simOpt.fracVert = 0.5;     % fraction of vertices initial set
simOpt.fracInpVert = 0.5;  % fraction of vertices input set
simOpt.inpChanges = 6;     % changes of input over time horizon

% random simulation
simRes = simulateRandom(tank,params,simOpt);


% Visualization -----------------------------------------------------------

dims = {[1,2],[3,4],[5,6]};

% plot different projections
for i = 1:length(dims)
    
    figure; hold on; box on;
    projDims = dims{i};

    % plot reachable sets
    hanParam = plot(RcontParam,projDims,'FaceColor',[.7 .7 .7],'EdgeColor','none');
    hanNoParam = plot(RcontNoParam,projDims,'FaceColor',[.5 .5 .5],'EdgeColor','none');
    
    % plot initial set
    plot(params.R0,projDims,'w','Filled',true,'EdgeColor','k');
  
    % plot simulation results
    plot(simRes,projDims);

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
    legend([hanParam,hanNoParam],'parametric','non-parametric');
end

completed = 1;

%------------- END OF CODE --------------
