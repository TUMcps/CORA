function res = testLong_nonlinParamSys_reach_01_tank
% testLong_nonlinParamSys_reach_01_tank - unit_test_function of
%    nonlinear reachability analysis with uncertain parameters
%
% Checks the solution of the nonlinearSys class for the 6 tank example with 
% uncertain parameters;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:
%    res = testLong_nonlinParamSys_reach_01_tank
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       30-June-2009
% Last update:   16-June-2011
%                16-August-2016
%                23-April-2020 (restructure params/options) 
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameters --------------------------------------------------------------

dim_x=6;
params.tFinal=400; %final time
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
    [2.2951563475236298; 2.1634660130913415; 2.0001713162063837; 1.8382964545932126; 1.6357089163639560; 1.3070344587514788], ...
    [3.8928903890554940; 3.7300482565925814; 3.4909185224659325; 3.2377181952313689; 3.1029347886879859; 3.2768262662390524]);

%final result
res = isequal(IH,IH_saved,1e-8);

% ------------------------------ END OF CODE ------------------------------
