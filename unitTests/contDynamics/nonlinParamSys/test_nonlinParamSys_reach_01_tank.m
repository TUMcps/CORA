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
params.tFinal=10; %final time
params.R0=zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim_x)]);
params.U=zonotope([0,0.005]);

%parameter intervals for reachability analysis
params.paramInt=interval(0.0148,0.015);


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

tankParam = nonlinParamSys(@tank6paramEq);


% Reachability Analysis ---------------------------------------------------
        
R = reach(tankParam,params,options);


% Numerical Evaluation ----------------------------------------------------
IH = interval(R.timeInterval.set{end});

%saved result
IH_saved = interval( ...
    [1.793326374412643, 3.425582765937606, 3.764452764459693, 2.115897469651695, 8.712849426635310, 4.398524309539285]', ...
    [2.262073865676699, 3.903416158189166, 4.183332823383807, 2.542076505174605, 9.273530229545282, 4.853518289475275]');

%final result
res = isequal(IH,IH_saved,1e-8);

% ------------------------------ END OF CODE ------------------------------
