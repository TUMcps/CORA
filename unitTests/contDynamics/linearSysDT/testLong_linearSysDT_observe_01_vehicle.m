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
% load vehicleModel_dim2 vehicle params options
[vehicle,params,options] = aux_load_model();

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
    
    % obtain enclosing intervals ---
    % VolMin
    if strcmp(options.alg,'VolMin-A') 
        IH_saved = [];
    elseif strcmp(options.alg,'VolMin-B') 
        IH_saved = interval( ...
            [-0.3707036241559320; 0.2870340937441627], ...
            [0.3897576159793223; 0.3003674270774960]);
        % FRad
    elseif strcmp(options.alg,'FRad-A') 
        IH_saved = interval( ...
            [-0.3090888964793908; 0.2849144661504785], ...
            [0.3031670983142311; 0.3022650442095819]);
    elseif strcmp(options.alg,'FRad-B') 
        IH_saved = interval( ...
            [-0.3118833062534078; 0.2859011229268012], ...
            [0.3052996942524376; 0.3013822318790924]);
    elseif strcmp(options.alg,'FRad-C') 
        IH_saved = interval( ...
            [-0.3118742761991223; 0.1358019644769430], ...
            [0.3048920303089872; 0.4357647602754841]);
        % PRad
    elseif strcmp(options.alg,'PRad-A') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-B') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-C') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-D') 
        IH_saved = [];
    elseif strcmp(options.alg,'PRad-E') 
        IH_saved = [];
        % Other
    elseif strcmp(options.alg,'Nom-G') 
        IH_saved = [];
    elseif strcmp(options.alg,'Hinf-G') 
        IH_saved = [];
        % ESO
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
        % CZN
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
    assert(isequal(IH,IH_saved,accuracy));
end

% final result
res = true;

end


% Auxiliary functions -----------------------------------------------------

function [vehicle,params,options] = aux_load_model()

% generated semi-automatically using printSystem, printSet, etc. 
% from saved .mat file

% build system
name = 'vehicle';
A = [ 0.9922748449450575 -0.0080030620640856 ; 0.2889726943723335 0.9299044703341978 ];
B = [ 0.0222329772960317 ; 0.2473223577905438 ];
c = [ 0.0000000000000000 ; 0.0000000000000000 ];
C = [ 0.0000000000000000 15.0000000000000000 ];
D = [];
k = [];
E = zeros(size(A,1),2);
F = zeros(size(C,1),1);
dt = 0.0100000000000000;
vehicle = linearSysDT(name,A,B,c,C,D,k,E,F,dt);

% build params
params = struct( ...
    'tStart', 0, ...
    'tFinal', 1, ...
    'R0', zonotope([ 0.0000000000000000 ; 0.0000000000000000 ], [ 1.0000000000000000 0.0000000000000000 ; 0.0000000000000000 1.0000000000000000 ]), ...
    'V', zonotope(0, 0.1000000000000000), ...
    'W', zonotope([ 0.0000000000000000 ; 0.0000000000000000 ], [ 0.0050000000000000 0.0000000000000000 ; 0.0000000000000000 0.0550000000000000 ]), ...
    'u', [ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000, 0.1000000000000000 ], ...
    'y', [ -4.7375685463570560, -3.6676088750789422, -3.0419374950055382, -2.6933739883929393, -2.4221039194812128, -2.3798534352542227, -1.8555146455315428, -2.1932338589400682, -1.7808281356382221, -1.4473209813924710, -1.4002972344402584, -1.1455625616985190, -0.5500999786962857, -0.4513452424928416, -0.0411425774708498, 0.4418651384431557, 0.5784758414529456, 0.2195897911424242, 0.5401026419837902, 0.6277260289661587, 0.2060559148172058, 0.2501216365593224, 0.5412751626338884, 0.3345866195444345, -0.0163675718377149, 0.3041195375290853, 0.6793302942259298, 1.0164514627077834, 1.5738742958676213, 1.1724695324159828, 1.2437964754547945, 1.4010164564118046, 1.5795757981189424, 1.6822370528578316, 1.4907111533152526, 1.5912651575329291, 1.6676719901865342, 1.6246405205105545, 1.9206853757952200, 1.4826667073424762, 1.1458610234447988, 1.2147916304264652, 1.2459390042675154, 0.8630460052549126, 0.6749731507225729, 0.9315718651089427, 1.4002607874103736, 1.8432981569635658, 1.7857096709612801, 2.0842649156005839, 1.9716635665952942, 2.0465576459276948, 2.3817969588476116, 2.7630618512163734, 3.1342017849807822, 3.3016675113514089, 3.1422060363039130, 3.0701264740032967, 3.2212536805204972, 3.1910453471250406, 3.6264482658971740, 4.0673519507826787, 3.5609837025049207, 3.7425418567197917, 4.2403880551103201, 3.8270923382285673, 4.2128693562095618, 4.9774380197745201, 4.9505775642274701, 5.0554870593120622, 4.4569108852355868, 4.5707653156727224, 4.8158225602269873, 4.4462950480472347, 4.8582327274512078, 4.6474673663467936, 4.8687243928277431, 4.3588351531233824, 4.5323646677790874, 4.3581118320076850, 4.5251743223319005, 4.1830309488322577, 4.5909334272268216, 4.3051590526590822, 4.8429703169499403, 4.8598948095052599, 4.6303192539380609, 4.6541280740345670, 4.7733859613167429, 4.5208540210281258, 4.6680710054598302, 4.3752270867401881, 4.3467381617989442, 4.4074930088629145, 4.0852651327984137, 3.9391444246783713, 3.9147287045366359, 3.9075387392212835, 4.2306228022559260, 4.4055114061624394 ] ...
);

% build options
options = struct( ...
    'timeStep', 0.0100000000000000, ...
    'zonotopeOrder', 100, ...
    'reductionTechnique', 'pca', ...
    'points', 1, ...
    'p_conf', 0.9990000000000000 ...
);

end

% ------------------------------ END OF CODE ------------------------------
