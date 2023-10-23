function res = test_linearSysDT_observe_02_2dim
% test_linearSysDT_observe_02_2dim - unit test function for guaranteed
% state estimation of linear discrete-time systems.
%
% Checks the solution of the linearSysDT class for a two-dimensional example 
% from Sec. 7.1 of [1]; It is checked whether the enclosing interval of the 
% final observed set is close to an interval provided by a previous 
% solution that has been saved
%
% Syntax:
%    res = test_linearSysDT_observe_02_2dim()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed
%        state estimation by zonotopes. Automatica, 41(6):1035–1043, 2005.

% Authors:       Matthias Althoff
% Written:       19-November-2020
% Last update:   04-June-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Load vehicle model
load twoDimSys sys params options

% save original parameters
paramsOrig = params;

% Set of evaluated estimators
Estimator = {
    'VolMin-B' 
    'FRad-A' 
    'FRad-B' 
    %'PRad-A' % tested in testSDPT3_linearSysDT_observe_02_2dim()
    %'PRad-B' % tested in testSDPT3_linearSysDT_observe_02_2dim()
    %'PRad-C' % tested in testSDPT3_linearSysDT_observe_02_2dim()
    %'PRad-D' % tested in testSDPT3_linearSysDT_observe_02_2dim()
    'CZN-A'
    'CZN-B'
    'ESO-A'
    'ESO-B'
    'FRad-C' 
    %'PRad-E' % tested in testSDPT3_linearSysDT_observe_02_2dim()
    %'Nom-G' % tested in testSDPT3_linearSysDT_observe_02_2dim()
    %'ESO-C' % tested in testSDPT3_linearSysDT_observe_02_2dim()
    %'ESO-D' % tested in testSDPT3_linearSysDT_observe_02_2dim()
    %'Hinf-G' % tested in testSDPT3_linearSysDT_observe_02_2dim()
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
    estSet = observe(sys,params,options);
    
    % enclose last estimated set by interval
    IH = interval(estSet.timePoint.set{end});
    
    % obtain enclosing intervals
    if strcmp(options.alg,'VolMin-A') 
        IH_saved = [];
    elseif strcmp(options.alg,'VolMin-B') 
        IH_saved = interval( ...
            [-0.1148743303581211; -0.3792153944480451], ...
            [0.1411256696418790; 0.2127846055519550]);
    elseif strcmp(options.alg,'FRad-A') 
        IH_saved = interval( ...
            [-0.1252836718233518; -0.3590605415736482], ...
            [0.1755077905850537; 0.3031460898598283]);
    elseif strcmp(options.alg,'FRad-B') 
        IH_saved = interval( ...
            [-0.1366939546269597; -0.3957333713260606], ...
            [0.1887846470750793; 0.3380097363827392]);
    elseif strcmp(options.alg,'PRad-A') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-B') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-C') 
        IH_saved = [];
    elseif strcmp(options.alg,'FRad-C') 
        IH_saved = interval( ...
            [-0.2942866270722150; -0.4084412460702273], ...
            [0.2767909177420897; 0.4113280380282030]);
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
            [-0.3846268741740866; -0.4863796623829761], ...
            [0.3847350367311811; 0.4863161984311742]);
    elseif strcmp(options.alg,'ESO-B')
        IH_saved = interval( ...
            [-3.7161600669997443; -4.6948479666014427], ...
            [3.7208223242092866; 4.6933762273831414]);
    elseif strcmp(options.alg,'ESO-C')
        IH_saved = [];   
    elseif strcmp(options.alg,'ESO-D')
        IH_saved = [];
    elseif strcmp(options.alg,'CZN-A')
        IH_saved = interval( ...
            [-0.1062309346120359; -0.2540037859300612], ...
            [0.1411256696418789; 0.2127846055519550]);
    elseif strcmp(options.alg,'CZN-B')
        IH_saved = interval( ...
            [-0.1062309346120359; -0.2540037859300612], ...
            [0.1411256696418789; 0.2127846055519550]);
    end

    %check if slightly bloated versions enclose each other
    resPartial(end+1) = isequal(IH,IH_saved,accuracy);
end

% final result
res = all(resPartial);
        
% ------------------------------ END OF CODE ------------------------------
