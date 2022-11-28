function res = example_nonlinearSysDT_reach_adaptive
% example_nonlinearSysDT_reach_adaptive - numerous examples for nonlinear
%    reachability analysis to test the adaptive parameter tuning approach
%
% Syntax:  
%    example_nonlinearSysDT_reach_adaptive
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
% 
% References:
%   ---

% Author:        Mark Wetzlinger
% Written:       21-September-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

res = true;

% System ------------------------------------------------------------------

n = 6; m = 1;
dt = 0.05;
sys = nonlinearSysDT(@tank6EqDT_V1,dt,n,m);


% Model parameters --------------------------------------------------------

params.tFinal = 40;
params.R0 = zonotope([2; 4; 4; 2; 10; 4],0.2*eye(n));
params.U = zonotope([0,0.005]);


% Reachability settings ---------------------------------------------------

options.alg = 'lin-adaptive';


% Reachability analysis ---------------------------------------------------

R = reach(sys,params,options);


% Visualization -----------------------------------------------------------

figure; hold on;
plot(R);


%------------- END OF CODE --------------
