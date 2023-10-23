function completed = example_linearSysDT_observe_01_vehicle()
% example_linearSysDT_observe_01_vehicle - example for guaranteed state
%    estimation of linear discrete-time systems; a corresponding unit test
%    function exists; shows the solution of the linearSysDT class for a
%    vehicle example.
%
% Syntax:
%    completed = example_linearSysDT_observe_01_vehicle()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false

% Authors:       Matthias Althoff
% Written:       30-April-2021
% Last update:   13-July-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Load vehicle model
load vehicleModel_dim6 vehicle params options simRes

% save original parameters
paramsOrig = params;

% Set of evaluated estimators
Estimator = {
    'FRad-A' 
    'FRad-B' 
    'ESO-A'
    'FRad-C'  
    'Nom-G' 
    'ESO-D'
    'Hinf-G' 
    };

% set solver
options.solver = 'sdpt3';
% ensure that solver is on path
if ~isSolverInstalled('sdpt3')
    throw(CORAerror('CORA:noSuitableSolver','sdpt3'));
end

%% perform evaluation
% loop over estimators
for iEst = 1:length(Estimator)

    % set algorithm
    estName = Estimator{iEst};
    options.alg = estName;
    
    % reset parameters
    params = paramsOrig;

    %% Initial sets, disturbance sets, and noise sets
    % if ellipsoids are required
    if any(strcmp(estName,{'ESO-A','ESO-B','ESO-C','ESO-D','ESO-E'}))
        params.R0 = ellipsoid(params.R0,'i:norm'); % inscribe ellipsoid
        params.W = ellipsoid(params.W,'i:norm'); % inscribe ellipsoid
        params.V = ellipsoid(params.V,'i:norm'); % inscribe ellipsoid
    % if constrained zonotopes are required
    elseif any(strcmp(estName,{'CZN-A','CZN-B'}))
        params.R0 = conZonotope(params.R0);
    end

    % perform set-based estimation 
    estSet = observe(vehicle,params,options);
    %estSetAll = evaluateSingleObserver(vehicle,params,options); % <-- uncomment for detailed evaluation
    %estSet = estSetAll.EstStates;
    
    % plot results
    for iDim = 1:2
        figure; hold on;
        % plot time elapse
        plotOverTime(estSet,iDim);
        % plot simulation
        plotOverTime(simRes,iDim);

        % label plot
        xlabel('t');
        ylabel(['x_{',num2str(iDim),'}']);
    end
end

% example completed
completed = true;


% ------------------------------ END OF CODE ------------------------------
