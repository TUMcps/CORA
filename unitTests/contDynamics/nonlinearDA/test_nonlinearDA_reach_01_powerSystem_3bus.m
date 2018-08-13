function res = test_nonlinearDA_reach_01_powerSystem_3bus()
% test_nonlinearDA_reach_01_powerSystem_3bus - unit_test_function of 
% nonlinear-differntial-algebraic reachability analysis
%
% Checks the solution of the nonlinearDASys class for a 3 bus power system example;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:  
%    res = test_nonlinearDA_reach_01_powerSystem_3bus()
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
% Written:      04-August-2016
% Last update:  19-August-2016
% Last revision:---


%------------- BEGIN CODE --------------

options.tensorOrder = 1;

%specify continuous dynamics-----------------------------------------------
powerDyn = nonlinDASys(2,6,2,@bus3Dyn,@bus3Con,options); %initialize power system dynamics
%--------------------------------------------------------------------------

%set options --------------------------------------------------------------
options.tStart = 0; %start time
options.tFinal = 1; %final time
options.x0 = [380; 0.7]; %initial state for simulation (consistemcy of algebraic state is taken care of automatically)
options.y0guess = [ones(0.5*powerDyn.nrOfConstraints, 1); zeros(0.5*powerDyn.nrOfConstraints, 1)];
options.R0 = zonotope([options.x0,diag([0.1, 0.01])]); %initial state for reachability analysis
options.uTrans = [1; 0.4];
options.U = zonotope([zeros(2,1),diag([0, 0.1*options.uTrans(2)])]);

%options.timeStep=0.01; %time step size for reachable set computation
options.timeStep=0.05; %time step size for reachable set computation
options.taylorTerms=6; %number of taylor terms for reachable sets
options.zonotopeOrder=200; %zonotope order
options.errorOrder=1.5;
options.polytopeOrder=2; %polytope order
options.reductionTechnique='girard';

options.originContained = 0;
options.reductionInterval = 1e5;
options.advancedLinErrorComp = 0;

options.maxError = [0.5; 0];
options.maxError_x = options.maxError;
options.maxError_y = 0.005*[1; 1; 1; 1; 1; 1];
%--------------------------------------------------------------------------

%compute reachable set 
Rcont = reach(powerDyn, options);

IH = interval(Rcont{end}{1});

%saved result
IH_saved = interval( ...
    [379.8304886258488864; 0.7548200082864], ...
    [381.7849447200206328; 0.7976176106355]);
        
%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_saved,1+1e-8));
res_2 = (IH_saved <= enlarge(IH,1+1e-8));

%final result
res = res_1*res_2;


%------------- END OF CODE --------------
        