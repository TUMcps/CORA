function completed = example_linearSysDT_observe_03_tank
% example_linearSysDT_observe_03_tank - example for guaranteed state
%    estimation of linear discrete-time systems from a unit test; shows the
%    solution of the linearSysDT class for a linearized tank example from [1].
%
% Syntax:
%    completed = example_linearSysDT_observe_03_tank
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% Reference:
%    [1] M. Althoff. Guaranteed State Estimation in CORA 2021. 
%        Proc. of the 8th International Workshop on Applied Verification 
%        for Continuous and Hybrid Systems, 2021.

% Authors:       Matthias Althoff
% Written:       12-June-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% load pedestrian model
load tankModel_lin_dim30 tank params options simRes

% observe
options.alg = 'FRad-C'; % observer approach
estSet = observe(tank,params,options);

% plot results
for iDim = 1:6
    figure; hold on;
    % plot time elapse
    plotOverTime(estSet,iDim);
    % plot simulation
    plotOverTime(simRes,iDim);

    % label plot
    title(options.alg);
    xlabel('t');
    ylabel(['x_{',num2str(iDim),'}']);
end

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
