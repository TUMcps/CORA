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
    [1.8058045990018012; 3.6440983425832627; 3.7944713775483399; 1.9523957411464474; 9.3414755040297610; 4.093338960605131], ...
    [2.2288260674127578; 4.0564919585455304; 4.1956261245334243; 2.3447014709840017; 9.7625790921187576; 4.485806366662812]);

IH_ti_true = interval( ...
    [1.7699897800185731; 3.6289355172711657; 3.7809745599505105; 1.7855046102852854; 9.3283653437240535; 3.790560381830072], ...
    [2.2489708336491194; 4.2199053537722717; 4.2152858794344752; 2.3647958718060287; 10.2195740971783469; 4.503745894115245]);

linErrors_true = 1e-3*[0.206844016449741; 0.121882639568035; 0.060496104052267; 0.251662494792440; 0.245782253631191; 0.097124423301666];
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
