function completed = example_nonlinearSysDT_observe_01_tank()
% example_nonlinearSysDT_observe_01_tank - example for guaranteed
% state estimation of nonlinear discrete-time systems from a unit test.
%
% Shows the solution of the nonlinearSysDT class for a tank example
% from [1].
%
% Syntax:  
%    completed = example_nonlinearSysDT_observe_01_tank
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
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] M. Althoff "Guaranteed State Estimation in CORA 2021", ARCH 2021
%
% Example: 
%    -
 
% Author:       Matthias Althoff
% Written:      14-Jun-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%% Load tank model
load tankModel_nonlin_dim30 tank params options simRes

% Set of evaluated estimators
Estimator = {
    'FRad-A' 
    'FRad-B' 
    'FRad-C' 
    };


%% perform evaluation
% loop over estimators
for iEst = 1:length(Estimator)

    % set algorithm
    estName = Estimator{iEst};
    options.alg = estName;
    
    % if constrained zonotopes or ellipsoids are used
    paramsNew = params;
    if any(strcmp(estName,{'ESO-A','ESO-B','ESO-C','ESO-D','ESO-E'}))
        paramsNew.R0 = ellipsoid(params.R0);
        paramsNew.W = ellipsoid(params.W);
        paramsNew.V = ellipsoid(params.V);
    elseif any(strcmp(estName,{'CZN-A','CZN-B'}))
        paramsNew.R0 = conZonotope(params.R0);
    end
    
    % evaluate observer
    estSet = observe(tank,paramsNew,options);

    % plot results
    if ~isempty(estSet)
        for dim = 1:6
            figure; hold on;
            % plot time elapse
            plotOverTime(estSet,dim,'FaceColor',[.6 .6 .6],'EdgeColor','none');
            % plot simulation
            plotOverTime(simRes,dim);

            % label plot
            title(estName);
            xlabel('t');
            ylabel(['x_{',num2str(dim),'}']);
        end
    end
end

% example completed
completed = 1;

%------------- END OF CODE --------------
