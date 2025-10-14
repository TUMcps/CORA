function res = test_nonlinDASys_simulate
% test_nonlinDASys_simulate - unit test for simulate
%
% Syntax:
%    res = test_nonlinDASys_simulate
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
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

% define sys --------------------------------------------------------------

states = 2;
inputs = 2;
f = @(x,y,u) x(1)+1+u(1);
g = @(x,y,u) (x(1)+1)*y(1) + 2;
sys = nonlinDASys(f,g);


% model parameters --------------------------------------------------------

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
[t,z,ind] = simulate(sys,params);
assert(size(t,2) == size(z,2))

% all checks ok
res = true;
end


% ------------------------------ END OF CODE ------------------------------
