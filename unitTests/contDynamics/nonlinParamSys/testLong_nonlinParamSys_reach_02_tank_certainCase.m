function res = testLong_nonlinParamSys_reach_02_tank_certainCase
% testLong_nonlinParamSys_reach_02_tank_certainCase - test of nonlinear
%    reachability analysis with uncertain parameters.
%    The difference compared to the uncertain case is that the parameter
%    value is fixed. Unlike the example in the class nonlinearSys, one can
%    change the parameter values using the options.paramInt.
%
% This example can be found in [1, Sec. 3.4.5] or in [2].
%
% Syntax:
%    res = testLong_nonlinParamSys_reach_02_tank_certainCase
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] M. Althoff, "Reachability analysis and its application to the
%        safety assessment of autonomous cars", Dissertation, TUM 2010.
%    [2] M. Althoff, O. Stursberg, and M. Buss. Reachability analysis 
%        of nonlinear systems with uncertain parameters using
%        conservative linearization. CDC 2008.

% Authors:       Matthias Althoff
% Written:       19-August-2016
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameters --------------------------------------------------------------

dim_x=6;
params.tFinal=400; %final time
params.R0=zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim_x)]);
params.U=zonotope([0,0.005]);


% Reachability Settings ---------------------------------------------------

options.timeStep=1;
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=10; %zonotope order
options.maxError = 1*ones(dim_x,1);
options.reductionInterval=1e3;
options.tensorOrder = 2;
options.alg = 'lin';


% System Dynamics with and without uncertain parameters -------------------

tankParam = nonlinParamSys(@tank6paramEq); %with uncertain parameters
tank = nonlinearSys(@tank6Eq); %without uncertain parameters


% Reachability Analysis ---------------------------------------------------

%compute reachable set of tank system with and without uncertain parameters
R_NoParam = reach(tank, params, options); %without uncertain parameters
params.paramInt=0.015; %parameter intervals for reachability analysis
options.intermediateTerms = 4;
R_Param = reach(tankParam,params,options); %with uncertain parameters


% Numerical Evaluation ----------------------------------------------------

%obtain interval hulls of both computation techniques
IHcontParam = interval(R_Param.timeInterval.set{end});
IHcontNoParam = interval(R_NoParam.timeInterval.set{end});

%final result
res = isequal(IHcontParam,IHcontNoParam,1e-8);

% ------------------------------ END OF CODE ------------------------------
