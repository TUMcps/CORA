function res = test_nonlinear_reach_01_tank(~)
% test_nonlinear_reach_01_tank - unit_test_function of nonlinear reachability analysis
%
% Checks the solution of the nonlinearSys class for the 6 tank example;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:  
%    res = test_nonlinear_reach_01_tank(~)
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
% Written:      21-July-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

dim=6;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=400; %final time
options.x0=[2; 4; 4; 2; 10; 4]; %initial state for simulation
options.R0=zonotope([options.x0,0.2*eye(dim)]); %initial state for reachability analysis

options.timeStep=4; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=50; %zonotope order
options.intermediateOrder=5;
options.reductionTechnique='girard';
options.errorOrder=1;
options.polytopeOrder=2; %polytope order
options.reductionInterval=1e3;
options.maxError = 1*ones(dim,1);

options.plotType='frame';
options.projectedDimensions=[1 2];

options.originContained = 0;
options.advancedLinErrorComp = 0;
options.tensorOrder = 2;
%--------------------------------------------------------------------------


%obtain uncertain inputs
options.uTrans = 0;
options.U = zonotope([0,0.005]); %input for reachability analysis

%specify continuous dynamics-----------------------------------------------
tank = nonlinearSys(6,1,@tank6Eq,options); %initialize tank system
%--------------------------------------------------------------------------


%compute reachable set using zonotopes
Rcont = reach(tank, options);

IH = interval(Rcont{end}{1});

%saved result
IH_saved = interval( ...
    [2.2939632391424079; 2.1709617663808629; 2.0160379729752376; 1.8476417858302909; 1.6416504727689167; 1.2970646473070846], ...
    [3.8374801442502351; 3.6739840853646228; 3.4336439419532145; 3.1915097908185457; 3.0610590827137658; 3.2471867050256562]);
        
%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_saved,1+1e-8));
res_2 = (IH_saved <= enlarge(IH,1+1e-8));

%final result
res = res_1*res_2;

%------------- END OF CODE --------------
