function res = test_nonlinearSysDT_observe_01_tank()
% test_nonlinearSysDT_observe_01_tank - unit_test_function for guaranteed
% state estimation of nonlinear discrete-time tank system.
%
% Checks the state estmation of the nonlinearSysDT class for a tank example 
% from [1]; It is checked whether the enclosing interval of the 
% final observed set is close to an interval provided by a previous 
% solution that has been saved
%
% Syntax:  
%    res = test_nonlinearSysDT_observe_01_tank
%
% Inputs:
%    -
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
% Written:      25-Mar-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%% Load tank model
load tankModel_nonlin_dim6 tank params options

% more robust case
params.tFinal = 10;
params.u = params.u(:,1:20);
params.y = params.y(:,1:20);

% Set of evaluated estimators
Estimator = {
    'FRad-A' 
    'FRad-B' 
    'FRad-C' 
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
    
    % enclose last estimated set by interval
    IH = interval(estSet.timePoint.set{end});
    
    % obtain enclosing intervals
    if strcmp(options.alg,'FRad-A')
        if params.tFinal == 10
            IH_saved = interval( ...
                [19.2976721922255834; 19.6139972122578357; 14.2446680848811198; 20.0892697254588484; 18.9581928439432339; 16.1726262416728765], ...
                [19.7876666310276192; 20.1035308207714820; 30.1688768049905462; 20.7224260892464685; 19.3788828254191650; 23.7131744879695638]);
        else
            IH_saved = interval( ...
                [17.6813648886554518; 18.2156906209731488; 21.2088894208727972; 22.3622857916818099; 22.4614208115655352; 19.2343954709766152], ...
                [18.3586147868548224; 18.9069663981664888; 23.8031380887066000; 23.2005173305999541; 23.1757768194797045; 22.9771198272602817]);
        end
    elseif strcmp(options.alg,'FRad-B')
        if params.tFinal == 10
            IH_saved = interval( ...
                [19.2712596426652674; 19.5913693494689234; 14.2056931305107454; 20.0980777032700715; 18.9587395711021074; 16.1529188853400001], ...
                [19.7900964406908493; 20.1029076589013762; 30.4543302458481833; 20.7268550247347747; 19.3789984819565717; 23.7268155612348224]);
        else
            IH_saved = interval( ...
                [17.6928498777203806; 18.2262865296234970; 21.2428263743206962; 22.3924573528288349; 22.4929791654838702; 19.2318792831541572], ...
                [18.3491195855406168; 18.8984370397519434; 23.7623779649680849; 23.1693183620486991; 23.1497672369610150; 22.9805998157471549]);
        end
    elseif strcmp(options.alg,'FRad-C')
        if params.tFinal == 10
            IH_saved = interval( ...
                [19.2882829662340818; 19.6049109089969669; 14.1375551339295349; 20.0683709752421144; 18.9552007160255158; 16.1715832824749270], ...
                [19.8034985374044119; 20.1189759828253436; 30.1814503670473009; 20.7408289024531278; 19.3772251039683709; 23.7132313182795755]);
        else
            IH_saved = interval( ...
                [17.6996609730186378; 18.2366127949558070; 21.3886795671867809; 22.3511769355149319; 22.5106975175051005; 19.2412346771689293], ...
                [18.3377371228686741; 18.8833803746827193; 23.6252998964590724; 23.2109899193938354; 23.1329874737624905; 22.9707010569532990]);
        end
    end

    %check if slightly bloated versions enclose each other
    resPartial(end+1) = (IH <= enlarge(IH_saved,1+accuracy));
    resPartial(end+1) = (IH_saved <= enlarge(IH,1+accuracy));
end

% final result
res = all(resPartial);

%------------- END OF CODE --------------
