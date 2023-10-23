function res = test_nonlinearSysDT_reach_02_tank
% test_nonlinearSysDT_reach_02_tank - unit_test_function of nonlinear
%    reachability analysis for discrete time systems. This unit test is the
%    discrete-time version of the unit test of the continuous-time tank 
%    system from [1]; It is checked whether the enclosing interval
%    of the final reachable set is close to an interval provided
%    by a previous solution that has been saved
%
% Syntax:
%    res = test_nonlinearSysDT_reach_02_tank
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
% Written:       25-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

dim_x=6;
params.tFinal = 100;
params.R0 = zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim_x)]);
params.U = zonotope([0,0.005]);


% Reachability Settings ---------------------------------------------------

% zonotope order: must be so high because the saved result was computed
%   at a time when the reduction-operation was not implemented
%   (max. order ~ 434)
options.zonotopeOrder = 450;
options.errorOrder = 1;

options.tensorOrder = 2;


% System Dynamics----------------------------------------------------------

% sampling time
dt = 0.5;

fun = @(x,u) tank6EqDT(x,u,dt);

tank = nonlinearSysDT('tankSystem',fun,dt); % initialize tank system


% Reachability Analysis ---------------------------------------------------

R = reach(tank,params,options);


% Numerical Evaluation ----------------------------------------------------

IH = interval(R.timePoint.set{end});

IH_saved = interval( ...
    [1.4029086368264769; 1.3837238402314385; 1.6839009122021129; 2.1098248033514708; 3.7063871491757219; 4.9884130326896310], ...
[1.9563699892722184; 1.9755297952057655; 2.2953017401028579; 2.6705455134974536; 4.3116242652586525; 5.4868053393744178]);

%final result
res = isequal(IH,IH_saved,1e-8);

% ------------------------------ END OF CODE ------------------------------
