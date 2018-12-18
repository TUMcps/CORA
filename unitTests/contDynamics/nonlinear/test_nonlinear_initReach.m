function res = test_nonlinear_initReach(~)
% test_nonlinear_initReach - unit_test_function for computing a single
% time intreval for nonlinear dynamics
%
% Checks initReach of the nonlinearSys class for the 6 tank example;
% It is checked whether partial reachable sets and the set of linearization
% errors are correctly obtained
%
% Syntax:  
%    res = test_nonlinear_initReach(~)
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      31-July-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

dim=6;

%set options --------------------------------------------------------------
options.uTrans = 0;
options.U = zonotope([0,0.005]); %input for reachability analysis
options.timeStep=4; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=50; %zonotope order
options.intermediateOrder=5;
options.reductionTechnique='girard';
options.maxError = 1*ones(dim,1);
options.originContained = 0;
options.advancedLinErrorComp = 0;
options.tensorOrder = 2;
%--------------------------------------------------------------------------

% obtain factors for reachability analysis
for i=1:(options.taylorTerms+1)
    %time step
    r = options.timeStep;
    %compute initial state factor
    options.factor(i)= r^(i)/factorial(i);    
end

%specify continuous dynamics-----------------------------------------------
tank = nonlinearSys(6,1,@tank6Eq,options); %initialize tank system
%--------------------------------------------------------------------------

%linearize system
R0=zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim)]); %initial state for reachability analysis
Rfirst = initReach(tank,R0,options);

% obtain interval hull of reachable set of first point in time
IH_tp = interval(Rfirst.tp{1}.set);
% obtain interval hull of reachable set of first time interval
IH_ti = interval(Rfirst.ti{1});
% obtain linearization errors
linErrors = Rfirst.tp{1}.error;

% provide ground truth-----------------------------------------------------
IH_tp_true = interval( ...
    [1.805794924492548; 3.643302925448510; 3.794026010097514; 1.951955268722969; 9.340994919972307; 4.092865541209832], ...
    [2.228835741922012; 4.057287375680283; 4.196071491984251; 2.345141943407481; 9.763059676176212; 4.486279786058111]);

IH_ti_true = interval( ...
    [1.769980105509319; 3.628140100136413; 3.780529192499684; 1.785064137861807; 9.327884759666599; 3.790086962434773], ...
    [2.248980508158373; 4.220700770907024; 4.215731246885301; 2.365236344229507; 10.220054681235801; 4.504219313510544]);

linErrors_true = 1e-3*[0.206863683556226; 0.314066832661960; 0.161658399976593; 0.353255589312750; 0.358487165091235; 0.209190685436450];
%--------------------------------------------------------------------------

%compare with obtained values
%check if slightly bloated versions enclose each other
res_tp_1 = (IH_tp <= enlarge(IH_tp_true,1+1e-8));
res_tp_2 = (IH_tp_true <= enlarge(IH_tp,1+1e-8));

res_ti_1 = (IH_ti <= enlarge(IH_ti_true,1+1e-8));
res_ti_2 = (IH_ti_true <= enlarge(IH_ti,1+1e-8));

res_error = (max(abs(linErrors - linErrors_true)) <= 1e-12);

%final result
res = res_tp_1*res_tp_2*res_ti_1*res_ti_2*res_error;

%------------- END OF CODE --------------
