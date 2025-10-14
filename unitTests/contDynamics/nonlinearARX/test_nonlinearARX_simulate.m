function res = test_nonlinearARX_simulate
% test_nonlinearARX_simulate - unit test for simulate
%
% Syntax:
%    res = test_nonlinearARX_simulate
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
f = @(y,u) [y(1,1) + u(1,1) - y(2,1); ...
              y(3,1) + u(2,1)*cos(y(1,1)); ...
              y(5,1) + u(4,1)*sin(y(1,1))];
dt = 0.25;
n_p = 2;
outputs = 3;
inputs = 2;
states = outputs * n_p;
sys = nonlinearARX(f,dt,outputs,inputs,n_p);

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
[t,~,~,y] = simulate(sys,params);
assert(size(t,2) == size(y,2))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
