function res = testLong_nonlinDASys_reach_01_powerSystem_3bus
% testLong_nonlinDASys_reach_01_powerSystem_3bus - unit_test_function of 
%    nonlinear-differntial-algebraic reachability analysis
%
% Checks the solution of the nonlinearDASys class for a 3 bus power system example;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:
%    res = testLong_nonlinDASys_reach_01_powerSystem_3bus
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       04-August-2016
% Last update:   19-August-2016
%                23-April-2020 (restructure params/options)
%                26-May-2022 (tighter computation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameters --------------------------------------------------------------

nrOfConstr = 6;
params.tFinal = 1; %final time
x0 = [380; 0.7]; %initial state for simulation (consistency of algebraic state is taken care of automatically)
params.y0guess = [ones(0.5*nrOfConstr,1); zeros(0.5*nrOfConstr,1)];
params.R0 = zonotope(x0,diag([0.1, 0.01])); %initial state for reachability analysis
params.U = zonotope([1; 0.4],diag([0, 0.04]));


% Reachability Settings ---------------------------------------------------

options.timeStep=0.05; %time step size for reachable set computation
options.taylorTerms=6; %number of taylor terms for reachable sets
options.zonotopeOrder=200; %zonotope order
options.tensorOrder = 2;
options.errorOrder = 3;
options.reductionInterval = 1e5;

options.maxError_x = [0.5; 0];
options.maxError_y = 0.005*[1; 1; 1; 1; 1; 1];


% System Dynamics ---------------------------------------------------------

powerDyn = nonlinDASys(@bus3Dyn,@bus3Con); %initialize power system dynamics


% Reachability Analysis ---------------------------------------------------

R = reach(powerDyn, params, options);


% Numerical Evaluation ----------------------------------------------------

IH = interval(R.timeInterval.set{end});

% saved result (computed using linError_mixed_noInt.m
% for linearization error computation)
IH_saved = interval( ...
    [379.8580735602910181; 0.7554235647684937], ...
    [381.7573146744116457; 0.7970152076943415]);
        
% compare results
res = isequal(IH,IH_saved,1e-8);

% ------------------------------ END OF CODE ------------------------------
