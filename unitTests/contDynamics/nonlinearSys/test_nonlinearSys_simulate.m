function res = test_nonlinearSys_simulate
% test_nonlinearSys_simulate - unit test for simulate
%
% Syntax:
%    res = test_nonlinearSys_simulate
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow
% Written:       28-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define system
fun = @(x,u) [x(2)+u(1); (1-x(1)^2)*x(2)-x(1)+u(2)];
sys = nonlinearSys('vanDerPol',fun);
inputs = 2;
states = 2;

% model parameters
% time horizon
params.tFinal = 1;
dt_steps = 10;
% initial set
params.R0 = zonotope(10*ones(states,1),0.5*diag(ones(states,1)));
% initial state
params.x0 = center(params.R0);
% input vector
u_m = randn(inputs,dt_steps+1);

% simulate ----------------------------------------------------------------
params.u = u_m;
[t,x] = simulate(sys,params);
assert(size(t,2) == size(x,2))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
