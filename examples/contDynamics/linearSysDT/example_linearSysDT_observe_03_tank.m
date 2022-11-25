function completed = example_linearSysDT_observe_03_tank
% example_linearSysDT_observe_03_tank - example for guaranteed
% state estimation of linear discrete-time systems from a unit test.
%
% Shows the solution of the linearSysDT class for a linearized tank example
% from [1].
%
% Syntax:  
%    completed = example_linearSysDT_observe_03_tank
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] M. Althoff. Guaranteed State Estimation in CORA 2021. 
%        Proc. of the 8th International Workshop on Applied Verification 
%        for Continuous and Hybrid Systems, 2021.
%
% Example: 
%    -
 
% Author:       Matthias Althoff
% Written:      12-Jun-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%% Load pedestrian model
load tankModel_lin_dim30 tank params options simRes

% observe
options.alg = 'FRad-C'; % observer approach
estSet = observe(tank,params,options);

% plot results
for dim = 1:6
    figure; hold on;
    % plot time elapse
    plotOverTime(estSet,dim,'FaceColor',[.6 .6 .6],'EdgeColor','none');
    % plot simulation
    plotOverTime(simRes,dim);

    % label plot
    title(options.alg);
    xlabel('t');
    ylabel(['x_{',num2str(dim),'}']);
end

% example completed
completed = 1;

%------------- END OF CODE --------------
