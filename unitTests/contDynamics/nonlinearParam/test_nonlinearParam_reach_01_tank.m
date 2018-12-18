function res = test_nonlinearParam_reach_01_tank()
% test_nonlinearParam_reach_01_tank - unit_test_function of nonlinear
% reachability analysis with uncertain parameters
%
% Checks the solution of the nonlinearSys class for the 6 tank example with 
% uncertain parameters;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:  
%    res = test_nonlinearParam_reach_01_tank()
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
% Written:      30-June-2009
% Last update:  16-June-2011
%               16-August-2016
% Last revision:---


%------------- BEGIN CODE --------------

dim=6;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=400; %final time
options.x0=[2; 4; 4; 2; 10; 4]; %initial state for simulation
options.R0=zonotope([options.x0,0.2*eye(dim)]); %initial state for reachability analysis
options.timeStep=1;
options.taylorTerms=4; %number of taylor terms for reachable sets
options.intermediateOrder = options.taylorTerms;
options.zonotopeOrder=10; %zonotope order
options.reductionTechnique='girard';
options.maxError = 1*ones(dim,1);
options.reductionInterval=1e3;
options.tensorOrder = 1;

options.advancedLinErrorComp = 0;

options.u=0; %input for simulation
options.U=zonotope([0,0.005]); %input for reachability analysis
options.uTrans=0; %has to be zero for nonlinear systems!!

options.p=0.015; %parameter values for simulation
options.paramInt=interval(0.0148,0.015); %parameter intervals for reachability analysis
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
tankParam = nonlinParamSys(6,1,1,@tank6paramEq,options); %initialize system
%tankParam = nonlinParamSys(6,1,1,@tank6paramEq,options.maxError,options); %initialize system
%--------------------------------------------------------------------------
        
%compute reachable set
Rcont = reach(tankParam,options);

%compute enclosing interval
IH = interval(Rcont{end}{1});

%saved result
IH_saved = interval( ...
    [2.2951563475236298; 2.1634660130913415; 2.0001713162063837; 1.8382964545932126; 1.6357089163639560; 1.3070344587514788], ...
    [3.8928903890554940; 3.7300482565925814; 3.4909185224659325; 3.2377181952313689; 3.1029347886879859; 3.2768262662390524]);

%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_saved,1+1e-8));
res_2 = (IH_saved <= enlarge(IH,1+1e-8));

%final result
res = res_1*res_2;


%------------- END OF CODE --------------