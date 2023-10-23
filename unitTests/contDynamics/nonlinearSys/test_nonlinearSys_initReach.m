function res = test_nonlinearSys_initReach
% test_nonlinearSys_initReach - unit_test_function for computing a single
%    time interval reachable set for nonlinear dynamics:
%    Checks initReach of the nonlinearSys class for the 6 tank example;
%    It is checked whether partial reachable sets and the set
%    of linearization errors are correctly obtained
%
% Syntax:
%    res = test_nonlinearSys_initReach
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       31-July-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Model Parameters --------------------------------------------------------

dim_x = 6;
% initial set for reachability analysis
params.R0=zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim_x)]);
% input for reachability analysis
params.U = zonotope([0,0.005]);
params.tFinal = 4;

% Reachability Settings ---------------------------------------------------

options.timeStep=4;
options.taylorTerms=4;
options.zonotopeOrder=50;
options.alg = 'lin';
options.tensorOrder = 2;


% System Dynamics ---------------------------------------------------------

tank = nonlinearSys(@tank6Eq);

% options check
options = validateOptions(tank,'reach',params,options);

% compute derivations (explicitly, since reach-function is not called)
derivatives(tank,options);

% obtain factors for reachability analysis
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end

% comupute only first step
Rfirst = initReach(tank,options.R0,options);

% obtain interval hull of reachable set of first point in time
IH_tp = interval(Rfirst.tp{1}.set);
% obtain interval hull of reachable set of first time interval
IH_ti = interval(Rfirst.ti{1});
% obtain linearization errors
linErrors = Rfirst.tp{1}.error;


% provide ground truth ----------------------------------------------------
IH_tp_true = interval( ...
    [1.805794924492548; 3.643302925448510; 3.794026010097514; 1.951955268722969; 9.340994919972307; 4.092865541209832], ...
    [2.228835741922012; 4.057287375680283; 4.196071491984251; 2.345141943407481; 9.763059676176212; 4.486279786058111]);

IH_ti_true = interval( ...
    [1.769980105509319; 3.628140100136413; 3.780529192499684; 1.785064137861807; 9.327884759666599; 3.790086962434773], ...
    [2.248980508158373; 4.220700770907024; 4.215731246885301; 2.365236344229507; 10.220054681235801; 4.504219313510544]);

linErrors_true = 1e-3*[0.206863683556226; 0.314066832661960; 0.161658399976593; 0.353255589312750; 0.358487165091235; 0.209190685436450];
% -------------------------------------------------------------------------

%final result
res = isequal(IH_tp,IH_tp_true,1e-8) && isequal(IH_ti,IH_ti_true,1e-8) ...
    && compareMatrices(linErrors,linErrors_true,1e-12);

% ------------------------------ END OF CODE ------------------------------
