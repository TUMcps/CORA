function res = testSDPT3_linearSysDT_observe_08_pedestrian()
% testSDPT3_linearSysDT_observe_08_pedestrian - unit_test_function for #
% guaranteed state estimation of linear discrete-time systems using the 
% SDPT3 solver (version 4.0).
%
% Checks the solution of the linearSysDT class for a pedestrian example;
% It is checked whether the enclosing interval of the final observed set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:
%    res = testSDPT3_linearSysDT_observe_08_pedestrian
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


%% Load pedestrian model
load pedestrianModel pedestrian params options


% Set of evaluated estimators
Estimator = {
    %'VolMin-B' % exclude for time saving of unit test
    %'FRad-A' % already tested in other unit test
    %'FRad-B' % already tested in other unit test
    'PRad-A' 
    'PRad-B' 
    'PRad-C' 
    'PRad-D' 
    %'CZN-A' % exclude for time saving of unit test
    %'CZN-B' % exclude for time saving of unit test
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
    estSet = observe(pedestrian,paramsNew,options);

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
            [1.0461829110338228; 2.1207241035823090; -1.0779763745894593; -0.4350782009753995], ...
            [1.3278380004635679; 2.4031365801691820; 1.2803330130003088; 1.9396781965441445]);
    elseif strcmp(options.alg,'PRad-B') 
        IH_saved = interval( ...
            [1.0370891373536133; 2.1110530212450813; -1.0360456003917879; -0.3980156282067054], ...
            [1.3367761434280636; 2.4117563165281322; 1.2365886092818661; 1.8925790574788914]);
    elseif strcmp(options.alg,'PRad-C') 
        IH_saved = interval( ...
            [1.0910512588820420; 2.1518210973083365; -1.1676141896315524; -0.6868646749564523], ...
            [1.3100607563815541; 2.3701215503011417; 1.4650347848631451; 2.0384371021192083]);
    elseif strcmp(options.alg,'FRad-C') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-D') 
        IH_saved = interval( ...
            [1.0993367738244246; 2.1661925142490950; -3.2126562053335554; -3.3849656968165145], ...
            [1.2999326462441796; 2.3668097391690472; 3.4833603242721614; 3.5537116116582665]);
    elseif strcmp(options.alg,'PRad-E') 
        IH_saved = interval( ...
            [1.0268347535869327; 2.1017257415132553; -1.0688025669859988; -0.4219264102185143], ...
            [1.3429774309843896; 2.4193540091223893; 1.2406407049675199; 1.9043020477561812]);
    elseif strcmp(options.alg,'Nom-G') 
        IH_saved = interval( ...
            [1.0825077416789488; 2.1024381408045127; -2.9502360499337312; -3.0131798798802896], ...
            [1.3453023181994352; 2.3676650157787260; 3.2681384852132807; 3.4457432052071404]);
    elseif strcmp(options.alg,'Hinf-G') 
        IH_saved = interval( ...
            [1.0468928642114645; 2.1321546573653563; -1.1738056318548948; -0.6691467301715885], ...
            [1.3069453761039502; 2.3944366253560241; 1.4503975578092114; 2.0521259320762888]);
    elseif strcmp(options.alg,'ESO-A')
        IH_saved = [];
    elseif strcmp(options.alg,'ESO-B')
        IH_saved = [];
    elseif strcmp(options.alg,'ESO-C')
        IH_saved = [];   
    elseif strcmp(options.alg,'ESO-D')
        IH_saved = interval( ...
            [-2.5125811063316421; -1.4547751412865555; -3.9944018412771549; -3.5340821447758728], ...
            [4.9059719923565881; 5.9637779572046163; 4.3177002187745597; 4.7780199151939708]);
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
