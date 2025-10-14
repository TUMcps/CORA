function res = test_linProbSys_simulate
% test_linProbSys_simulate - unit test for simulate
%
% Syntax:
%    res = test_linProbSys_simulate
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

% define linearSys --------------------------------------------------------

states = 2;
inputs = 2;
A = [-1 -4; 4 -1];
B = eye(2);
C = 0.7*eye(2);
sys = linProbSys(A,B,C);


% model parameters --------------------------------------------------------

% time horizon
params.tFinal = 1;
options.timeStep = 0.5;
% initial set
params.R0 = zonotope(10*ones(states,1),0.5*diag(ones(states,1)));
% initial state
params.x0 = center(params.R0);
% input vector (currently only constant input supported)
u_m = randn(inputs,1);

% simulate ----------------------------------------------------------------

params.u = u_m;
[t,x] = simulate(sys,params,options);
assert(size(t,2) == size(x,2))

% all checks ok
res = true;
end


% ------------------------------ END OF CODE ------------------------------
