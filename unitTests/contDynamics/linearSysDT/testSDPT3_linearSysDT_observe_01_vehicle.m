function res = testSDPT3_linearSysDT_observe_01_vehicle()
% testSDPT3_linearSysDT_observe_01_vehicle - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the 
% SDPT3 solver (version 4.0).
%
% Checks the solution of the linearSysDT class for a vehicle example;
% It is checked whether the enclosing interval of the final observed set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:
%    res = testSDPT3_linearSysDT_observe_01_vehicle()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 

% Authors:       Matthias Althoff
% Written:       07-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Load vehicle model
load vehicleModel_dim2 vehicle params options

% save original parameters
paramsOrig = params;

% Set of evaluated estimators
Estimator = {
    %'VolMin-B' % already tested in other unit test
    %'FRad-A' % already tested in other unit test
    %'FRad-B' % already tested in other unit test
    'PRad-A' 
    'PRad-B'
    'PRad-C' 
    'PRad-D' 
    %'CZN-A' % already tested in other unit test
    %'CZN-B' % already tested in other unit test
    %'ESO-A' % already tested in other unit test
    %'ESO-B' % already tested in other unit test
    %'FRad-C' % already tested in other unit test
    'PRad-E' 
    'Nom-G' 
    %'ESO-C' % does not succeed due badly conditioned matrices Q, P
    'ESO-D' 
    'Hinf-G' 
    };

% init
resPartial = [];

% set solver
options.solver = 'sdpt3';

%% perform evaluation
     
% loop over estimators
for iEst = 7:length(Estimator)

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
        IH_saved = [];
    elseif strcmp(options.alg,'FRad-A') 
        IH_saved = [];
    elseif strcmp(options.alg,'FRad-B') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-A') 
        IH_saved = interval( ...
            [-0.2172044379198976; 0.2866142489042661], ...
            [0.2111830817523198; 0.3007198979629934]);
    elseif strcmp(options.alg,'PRad-B') 
        IH_saved = interval( ...
            [-0.2172044379198976; 0.2866142489042661], ...
            [0.2111830817523198; 0.3007198979629934]);
    elseif strcmp(options.alg,'PRad-C') 
        IH_saved = interval( ...
            [-0.2486260163569280; 0.2848870172565884], ...
            [0.2233598172233524; 0.3021460167067332]);
    elseif strcmp(options.alg,'FRad-C') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-D') 
        IH_saved = interval( ...
            [-0.2421170756779888; 0.2870199690493845], ...
            [0.2164874543144344; 0.3003841936544813]);
    elseif strcmp(options.alg,'PRad-E') 
        IH_saved = interval( ...
            [-0.2291314405234857; 0.1604798435187392], ...
            [0.2058417740150857; 0.4066150669340608]);
    elseif strcmp(options.alg,'Nom-G') 
        IH_saved = interval( ...
            [-0.2537397902841674; 0.1522072571970656], ...
            [0.2219587794172041; 0.4110926389881205]);
    elseif strcmp(options.alg,'Hinf-G') 
        IH_saved = interval( ...
            [-0.2528566876515858; 0.1330898829152801], ...
            [0.2012227176653183; 0.3856656111217768]);
    elseif strcmp(options.alg,'ESO-A')
        IH_saved = [];
    elseif strcmp(options.alg,'ESO-B')
        IH_saved = [];
    elseif strcmp(options.alg,'ESO-C')
        IH_saved = [];   
    elseif strcmp(options.alg,'ESO-D')
        IH_saved = interval( ...
            [-1.9091774607938188; -1.3394047313881137], ...
            [1.8746418398906504; 1.9012670147926500]);
    elseif strcmp(options.alg,'CZN-A')
        IH_saved = [];
    elseif strcmp(options.alg,'CZN-B')
        IH_saved = [];
    end

    %check if slightly bloated versions enclose each other
    resPartial(end+1) = isequal(IH,IH_saved,1e-6);
end

% final result
res = all(resPartial);

% ------------------------------ END OF CODE ------------------------------
