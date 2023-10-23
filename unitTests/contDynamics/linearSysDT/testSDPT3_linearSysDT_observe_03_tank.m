function res = testSDPT3_linearSysDT_observe_03_tank()
% testSDPT3_linearSysDT_observe_03_tank - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the 
% SDPT3 solver (version 4.0).
%
% Checks the solution of the linearSysDT class for a tank example;
% It is checked whether the enclosing interval of the final observed set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:
%    res = testSDPT3_linearSysDT_observe_03_tank
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
load tankModel_lin_dim30 tank params options

% Set of evaluated estimators
Estimator = {
    %'VolMin-B' % exclude for time saving of unit test
    %'FRad-A' % already tested in other unit test
    %'FRad-B' % already tested in other unit test
    %'PRad-A' % exclude for time saving of unit test
    %'PRad-B' % exclude for time saving of unit test
    %'PRad-C' % exclude for time saving of unit test
    %'PRad-D' % exclude for time saving of unit test
    %'CZN-A' % exclude for time saving of unit test
    %'CZN-B' % exclude for time saving of unit test
    %'ESO-A' % already tested in other unit test
    %'ESO-B' % already tested in other unit test
    %'FRad-C' % already tested in other unit test
    %'PRad-E' % exclude for time saving of unit test
    'Nom-G' 
    %'ESO-C' % exclude for time saving of unit test
    %'ESO-D' % exclude for time saving of unit test
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
    
    paramsNew.U = zonotope(0); % <--- needs to be removed after updating checkOptions
    
    % evaluate observer
    estSet = observe(tank,paramsNew,options);
    
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
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-B') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-C') 
        IH_saved = [];
    elseif strcmp(options.alg,'FRad-C') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-D') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-E') 
        IH_saved = [];
    elseif strcmp(options.alg,'Nom-G') 
        IH_saved = interval( ...
            [25.0746399912694642; 26.2900911106279800; 8.5463471407233840; 4.7929462731059065; 1.3570965158383361; 0.9245609802752499; 1.7006136314942113; 0.6259160497144953; -1.3226080426961921; -0.2943868677169681; -1.0147311059334472; -0.0256423698820479; 0.5682207711008320; 1.4026420580074173; -0.2882666334892989; 2.7959133917470718; 3.2085636119976293; 2.5938033881634750; 3.2003693846634520; 3.6918763309861973; 2.4259484599967105; 2.1087672772854509; 1.5823018430897438; 0.2253056396716457; 0.8253869720454425; 1.1562057264225063; 1.1751343867185502; 1.7801864850564506; 1.6897711543184366; -0.4340077176522699], ...
            [25.4908057325857200; 26.7012827782155888; 12.4226315621522616; 5.2089736591781453; 1.7630718997728232; 1.4217686830014085; 2.1371262549766881; 1.0297837994859205; 2.7095534959862952; 0.1223155324659134; -0.6094746214857001; 0.5350017722603975; 1.0070650063933855; 1.8065077999160639; 3.7555120680184158; 3.2126624749553896; 3.6150353412484737; 3.1255769806211378; 3.6364898459302952; 4.1040769707369504; 2.8382281853011766; 2.5218339605327538; 1.9858275701496495; 0.7611971699117170; 1.2635988009098724; 1.5688220365465317; 1.5877654471510356; 2.1923650450214387; 2.1016385432101212; 3.2453441995767935]);
    elseif strcmp(options.alg,'Hinf-G') 
        IH_saved = interval( ...
            [25.1577152768565426; 26.4442193439533568; 8.9829441008665061; 4.8015811330788756; 1.3999806104755388; 1.0004102683044531; 1.7878064277745542; 0.7257144817444836; -1.6325897755445609; -0.2179877381373467; -1.2596475217000396; -0.0090196398146163; 0.5509798944603309; 1.5063680256762626; 0.2272435411147828; 2.9211857299282808; 3.0449324971855956; 2.6287409096040930; 3.1685016314244723; 3.7894074185397253; 2.5536673889289268; 2.1436026353315500; 1.4910160917548836; 0.2638570804729945; 0.8606288220311487; 1.2094204482572204; 1.2780219524943748; 1.6675582375766762; 1.8647940634604634; -0.4622347771465969], ...
            [25.5572961864667292; 26.8562003296650715; 12.7921930872872256; 5.2270539785355767; 1.8046553120336086; 1.4127443321600379; 2.1985913483938049; 1.1335759016477911; 2.1795927761323464; 0.2068017118548641; -0.8555198955177323; 0.4557594946051987; 0.9623881763983305; 1.9142306804462819; 4.0394956507107178; 3.3459778092002646; 3.4496066628583812; 3.0705599918177073; 3.5788444674475470; 4.1938155923536060; 2.9567334172368587; 2.5473666950067995; 1.8951246405808961; 0.7099490481813026; 1.2722793826087786; 1.6131951061329928; 1.6817875543169483; 2.0713218704819756; 2.2683007806350437; 3.2511657902145195]);
    elseif strcmp(options.alg,'ESO-A')
        IH_saved = [];
    elseif strcmp(options.alg,'ESO-B')
        IH_saved = [];
    elseif strcmp(options.alg,'ESO-C')
        IH_saved = [];   
    elseif strcmp(options.alg,'ESO-D')
        IH_saved = [];
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
