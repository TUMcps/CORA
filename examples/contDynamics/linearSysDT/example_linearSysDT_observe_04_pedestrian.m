function completed = example_linearSysDT_observe_04_pedestrian
% example_linearSysDT_observe_04_pedestrian - example for guaranteed
% state estimation of linear discrete-time systems from a unit test.
%
% Provides the solution of the linearSysDT class for a pedestrian example
%
% Syntax:  
%    completed = example_linearSysDT_observe_04_pedestrian
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%    -
 
% Author:       Matthias Althoff
% Written:      03-Mar-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%% Load pedestrian model
load pedestrianModel pedestrian params options simRes


% Set of evaluated estimators
Estimator = {
    'VolMin-B' 
    'FRad-A' 
    'FRad-B' 
    'PRad-A' 
    'PRad-B' 
    'PRad-C' 
    'PRad-D' 
    'ESO-A'
    'ESO-B'
    'FRad-C' 
    'PRad-E' 
    'Nom-G' 
    'ESO-D'
    'Hinf-G' 
    };

% set solver
options.solver = 'sdpt3';


%% perform evaluation
% loop over estimators
for iEst = 1:length(Estimator)

    % set algorithm
    estName = Estimator{iEst}
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
    estSet = observe(pedestrian,paramsNew,options);

    % plot results
    if ~isempty(estSet)
        for dim = 1:4
            figure; hold on;
            % plot time elapse
            plotOverTime(estSet,dim,'FaceColor',[.6 .6 .6],'EdgeColor','none');
            % plot simulation
            plotOverTime(simRes,dim);

            % label plot
            xlabel('t');
            ylabel(['x_{',num2str(dim),'}']);
        end
    end
end

% example completed
completed = 1;

%------------- END OF CODE --------------
