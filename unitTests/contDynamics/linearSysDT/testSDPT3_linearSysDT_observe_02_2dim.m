function res = testSDPT3_linearSysDT_observe_02_2dim
% testSDPT3_linearSysDT_observe_02_2dim - unit test function for 
% guaranteed state estimation of linear discrete-time systems using the 
% SDPT3 solver (version 4.0).
%
% Checks the solution of the linearSysDT class for a two-dimensional example 
% from Sec. 7.1 of [1]; It is checked whether the enclosing interval of the 
% final observed set is close to an interval provided by a previous 
% solution that has been saved
%
% Syntax:
%    res = testSDPT3_linearSysDT_observe_02_2dim()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed
%        state estimation by zonotopes. Automatica, 41(6):1035–1043, 2005.

% Authors:       Matthias Althoff
% Written:       07-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Load vehicle model
load twoDimSys sys params options

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
    params.U = zonotope(0); % <--- needs to be removed after updating checkOptions
    
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
        IH_saved = [];
    elseif strcmp(options.alg,'FRad-A') 
        IH_saved = [];
    elseif strcmp(options.alg,'FRad-B') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-A') 
        IH_saved = interval( ...
            [-0.1641491575720032; -0.4022035529669543], ...
            [0.2245099106298341; 0.3341163184271425]);
    elseif strcmp(options.alg,'PRad-B') 
        IH_saved = interval( ...
            [-0.1641491575720032; -0.4022035529669543], ...
            [0.2245099106298341; 0.3341163184271425]);
    elseif strcmp(options.alg,'PRad-C') 
        IH_saved = interval( ...
            [-0.1145617857673594; -0.3101234942629514], ...
            [0.1604095301851725; 0.2512275601704513]);
    elseif strcmp(options.alg,'FRad-C') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-D') 
        IH_saved = interval( ...
            [-0.1150749893102779; -0.3731416616256365], ...
            [0.1429357144726702; 0.2147702994613390]);
    elseif strcmp(options.alg,'PRad-E') 
        IH_saved = interval( ...
            [-0.3179286260499298; -0.4757940345698530], ...
            [0.2925684431937721; 0.4814143751781071]);
    elseif strcmp(options.alg,'Nom-G') 
        IH_saved = interval( ...
            [-0.3006414281693829; -0.3971585808753912], ...
            [0.2638590955900897; 0.4176788743867478]);
    elseif strcmp(options.alg,'Hinf-G') 
        IH_saved = interval( ...
            [-0.2426101873789709; -0.4082011824647768], ...
            [0.3058912729587950; 0.3575579923788003]);
    elseif strcmp(options.alg,'ESO-A')
        IH_saved = [];
    elseif strcmp(options.alg,'ESO-B')
        IH_saved = [];
    elseif strcmp(options.alg,'ESO-C')
        IH_saved = [];   
    elseif strcmp(options.alg,'ESO-D')
        IH_saved = interval( ...
            [-1.6532157758372144; -2.4633737374923004], ...
            [1.6117166192534249; 2.4912774745301847]);
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
