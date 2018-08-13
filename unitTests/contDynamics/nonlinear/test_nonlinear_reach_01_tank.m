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
    [2.6166933508340877; 2.5277851118365846; 2.3715839274694237; 2.2060200588589631; 2.0636082283502004; 1.979956874998989], ...
    [3.5147500325585552; 3.3171607399089011; 3.0780979874590284; 2.8331315177898735; 2.6391013271324821; 2.564294477333752]);
        
%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_saved,1+1e-8));
res_2 = (IH_saved <= enlarge(IH,1+1e-8));

%final result
res = res_1*res_2;

%------------- END OF CODE --------------
