function res = test_nonlinParamSys_reach_01_tank
% test_nonlinParamSys_reach_01_tank - unit_test_function of nonlinear
%    reachability analysis with uncertain parameters
%
% Checks the solution of the nonlinearSys class for the 6 tank example with 
% uncertain parameters;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:  
%    res = test_nonlinParamSys_reach_01_tank
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%    -

% Author:       Matthias Althoff
% Written:      30-June-2009
% Last update:  16-June-2011
%               16-August-2016
%               23-April-2020 (restructure params/options) 
% Last revision:---

%------------- BEGIN CODE --------------


% Parameters --------------------------------------------------------------

dim_x=6;
params.tFinal=100; %final time
params.R0=zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim_x)]);
params.U=zonotope([0,0.005]);

params.paramInt=interval(0.0148,0.015); %parameter intervals for reachability analysis


% Reachability Settings ---------------------------------------------------

options.timeStep=1;
options.taylorTerms=4; %number of taylor terms for reachable sets
options.intermediateTerms = 4;
options.zonotopeOrder=10; %zonotope order
options.maxError = 1*ones(dim_x,1);
options.reductionInterval=1e3;
options.tensorOrder = 2;
options.alg = 'lin';

% System Dynamics ---------------------------------------------------------
tankParam = nonlinParamSys(@tank6paramEq); %initialize system


% Reachability Analysis ---------------------------------------------------
        
R = reach(tankParam,params,options);


% Numerical Evaluation ----------------------------------------------------
IH = interval(R.timeInterval.set{end});

%saved result
IH_saved = interval( ...
    [1.3349391829238719, 1.6195288465263205, 2.3126007171678191, 2.7486664833309211, 4.1540823219929308, 5.1469704676249135]', ...
    [2.0055927587397444, 2.3873588956536591,3.0825449137000458,3.3418651908699388,4.8034866551097544,5.6908886836716226]');
%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_saved,1+1e-8));
res_2 = (IH_saved <= enlarge(IH,1+1e-8));

%final result
res = res_1 && res_2;


%------------- END OF CODE --------------