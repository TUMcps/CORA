function res = testLong_linearSysDT_observe_01_vehicle()
% testLong_linearSysDT_observe_01_vehicle - unit_test_function for guaranteed
% state estimation of linear discrete-time systems.
%
% Checks the solution of the linearSysDT class for a vehicle example;
% It is checked whether the enclosing interval of the final observed set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:
%    res = testLong_linearSysDT_observe_01_vehicle()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Matthias Althoff
% Written:       08-September-2020
% Last update:   16-January-2021
%                04-June-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Load vehicle model
load vehicleModel_dim2 vehicle params options

% save original parameters
paramsOrig = params;

% Set of evaluated estimators
Estimator = {
    'VolMin-B' 
    'FRad-A' 
    'FRad-B' 
    %'PRad-A' % tested in testSDPT3_linearSysDT_observe_01_vehicle()
    %'PRad-B' % tested in testSDPT3_linearSysDT_observe_01_vehicle()
    %'PRad-C' % tested in testSDPT3_linearSysDT_observe_01_vehicle()
    %'PRad-D' % tested in testSDPT3_linearSysDT_observe_01_vehicle()
    'CZN-A'
    'CZN-B'
    'ESO-A'
    'ESO-B'
    'FRad-C' 
    %'PRad-E' % tested in testSDPT3_linearSysDT_observe_01_vehicle()
    %'Nom-G' % tested in testSDPT3_linearSysDT_observe_01_vehicle()
    %'ESO-C' % tested in testSDPT3_linearSysDT_observe_01_vehicle()
    %'ESO-D' % tested in testSDPT3_linearSysDT_observe_01_vehicle()
    %'Hinf-G' % tested in testSDPT3_linearSysDT_observe_01_vehicle()
    };

% init
resPartial = [];

% set accuracy
accuracy = 1e-6;

%% perform evaluation
% loop over estimators
for iEst = 1:length(Estimator)

    % set algorithm
    estName = Estimator{iEst};
    options.alg = estName;

    % reset parameters
    params = paramsOrig;
    
    % if ellipsoids are required
    if any(strcmp(estName,{'ESO-A','ESO-B','ESO-C','ESO-D','ESO-E'}))
        params.R0 = ellipsoid(params.R0,'inner:norm'); % inscribe ellipsoid
        params.W = ellipsoid(params.W,'inner:norm'); % inscribe ellipsoid
        params.V = ellipsoid(params.V,'inner:norm'); % inscribe ellipsoid
    % if constrained zonotopes are required
    elseif any(strcmp(estName,{'CZN-A','CZN-B'}))
        params.R0 = conZonotope(params.R0);
    end

    % run observer
    estSet = observe(vehicle,params,options);
    
    % enclose last estimated set by interval
    IH = interval(estSet.timePoint.set{end});
    
    % obtain enclosing intervals
    if strcmp(options.alg,'VolMin-A') 
        IH_saved = [];
    elseif strcmp(options.alg,'VolMin-B') 
        IH_saved = interval( ...
            [-0.3707036241559320; 0.2870340937441627], ...
            [0.3897576159793223; 0.3003674270774960]);
    elseif strcmp(options.alg,'FRad-A') 
        IH_saved = interval( ...
            [-0.3090888964793908; 0.2849144661504785], ...
            [0.3031670983142311; 0.3022650442095819]);
    elseif strcmp(options.alg,'FRad-B') 
        IH_saved = interval( ...
            [-0.3118833062534078; 0.2859011229268012], ...
            [0.3052996942524376; 0.3013822318790924]);
    elseif strcmp(options.alg,'PRad-A') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-B') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-C') 
        IH_saved = [];
    elseif strcmp(options.alg,'FRad-C') 
        IH_saved = interval( ...
            [-0.3118742761991223; 0.1358019644769430], ...
            [0.3048920303089872; 0.4357647602754841]);
    elseif strcmp(options.alg,'PRad-D') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-E') 
        IH_saved = [];
    elseif strcmp(options.alg,'Nom-G') 
        IH_saved = [];
    elseif strcmp(options.alg,'Hinf-G') 
        IH_saved = [];
    elseif strcmp(options.alg,'ESO-A')
        IH_saved = interval( ...
            [-0.2490582834001083; -0.7714450648968025], ...
            [0.2343518379626528; 1.4310874124628041]);
    elseif strcmp(options.alg,'ESO-B')
        IH_saved = interval( ...
            [-0.5455598011825439; -1.9773034896935415], ...
            [0.5293559999343277; 2.6167343485964136]);
    elseif strcmp(options.alg,'ESO-C')
        IH_saved = [];   
    elseif strcmp(options.alg,'ESO-D')
        IH_saved = [];
    elseif strcmp(options.alg,'CZN-A')
        IH_saved = interval( ...
            [-0.1808439583069218; 0.2870340937441625], ...
            [0.1636688236433201; 0.3003674270774961]);
    elseif strcmp(options.alg,'CZN-B')
        IH_saved = interval( ...
            [-0.1808439583069214; 0.2870340937441626], ...
            [0.1636688236433202; 0.3003674270774960]);
    end

    %check if slightly bloated versions enclose each other
    resPartial(end+1) = isequal(IH,IH_saved,accuracy);
end

% final result
res = all(resPartial);

% ------------------------------ END OF CODE ------------------------------
