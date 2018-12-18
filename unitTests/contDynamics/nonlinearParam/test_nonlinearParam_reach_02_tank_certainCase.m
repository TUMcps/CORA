function res = test_nonlinearParam_reach_02_tank_certainCase()
% test_nonlinearParam_reach_02_tank_certainCase - test of nonlinear
% reachability analysis with uncertain parameters. The difference compared to the uncertain case is
% that the parameter value is fixed. Unlike the example in the class
% nonlinearSys, one can change the parameter values using the
% options.paramInt.
%
% This example can be found in Sec. 3.4.5 of
% M. Althoff, “Reachability analysis and its application to the safety 
% assessment of autonomous cars”, Dissertation, Technische Universität 
% München, 2010, 
% http://nbn-resolving.de/urn/resolver.pl?urn:nbn:de:bvb:91-diss-20100715-963752-1-4
% 
% or in
%
% M. Althoff, O. Stursberg, and M. Buss. Reachability analysis of nonlinear 
% systems with uncertain parameters using conservative linearization. In 
% Proc. of the 47th IEEE Conference on Decision and Control, 
% pages 4042–4048, 2008
%
% Syntax:  
%    test_nonlinearParam_reach_02_tank_certainCase
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
% Written:      19-August-2016
% Last update:  ---
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
options.uTrans=0; 

options.p=0.015; %parameter values for simulation
options.paramInt=0.015; %parameter intervals for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics with and without uncertain parameters---------
tankParam = nonlinParamSys(6,1,1,@tank6paramEq,options); %with uncertain parameters
tank = nonlinearSys(6,1,@tank6Eq,options); %without uncertain parameters
%--------------------------------------------------------------------------
        
%compute reachable set of tank system with and without uncertain parameters
RcontParam = reach(tankParam,options); %with uncertain parameters
RcontNoParam = reach(tank, options); %without uncertain parameters

%obtain interval hulls of both computation techniques
IHcontParam = interval(RcontParam{end}{1});
IHcontNoParam = interval(RcontNoParam{end}{1});

%check if slightly bloated versions enclose each other
res_1 = (IHcontParam <= enlarge(IHcontNoParam,1+1e-8));
res_2 = (IHcontNoParam <= enlarge(IHcontParam,1+1e-8));

%final result
res = res_1*res_2;

%------------- END OF CODE --------------
