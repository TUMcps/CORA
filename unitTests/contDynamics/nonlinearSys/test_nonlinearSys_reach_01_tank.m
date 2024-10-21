function res = test_nonlinearSys_reach_01_tank
% test_nonlinearSys_reach_01_tank - unit_test_function of nonlinear
%    reachability analysis; Checks the solution of the nonlinearSys class
%    for the 6 tank example from [1]; It is checked whether the enclosing 
%    interval of the final reachable set is close to an interval provided
%    by a previous solution that has been saved
%
% Syntax:
%    res = test_nonlinearSys_reach_01_tank
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] M. Althoff, O. Stursberg, and M. Buss: Reachability Analysis of 
%        Nonlinear Systems with Uncertain Parameters using Conservative 
%        Linearization. Proc. of the 47th IEEE Conference on Decision and 
%        Control, 2008.

% Authors:       Matthias Althoff
% Written:       21-July-2016
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% Parameters --------------------------------------------------------------

dim_x=6;
params.tFinal=400; %final time
params.R0=zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim_x)]);
params.U = zonotope([0,0.005]);


% Reachability Settings ---------------------------------------------------

options.timeStep=4; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=50; %zonotope order

options.alg = 'lin';
options.tensorOrder = 2;


% System Dynamics----------------------------------------------------------

tank = nonlinearSys(@tank6Eq); %initialize tank system


% Reachability Analysis ---------------------------------------------------

R = reach(tank,params,options);


% Numerical Evaluation ----------------------------------------------------

IH = interval(R.timeInterval.set{end});

IH_saved = interval( ...
    [2.2936477632745289; 2.1706800635653218; 2.0156503777163191; 1.8474995639752518; 1.6412146827755991; 1.2961102589990721], ...
    [3.8377956209405877; 3.6742657882546377; 3.4340315367183956; 3.1916520121356604; 3.0614948721719735; 3.2481410926468546]);

%final result
assert(isequal(IH,IH_saved,1e-8));

% ------------------------------ END OF CODE ------------------------------
