function res = test_linearSys_simulate_time
% test_linearSys_simulate_time - unit test for simulate specifically
%    checking whether start time, end time, and steps are correct
%
% Syntax:
%    res = test_linearSys_simulate_time
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

% Authors:       Mark Wetzlinger
% Written:       05-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
n = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];

% output matrix
C = [1 1 0 0;
     0 -0.5 0.5 0];

% initialize different linearSys-objects
sys = linearSys(A,B,zeros(4,1),C);

% initial set
R0 = zonotope(10*ones(n,1),0.5*diag(ones(n,1)));
% initial state
params.x0 = center(R0);

% only time horizon
params.tFinal = 1;
[t,x] = simulate(sys,params);

% check if start and end time correct
res(end+1,1) = t(1) == 0 && withinTol(t(end),params.tFinal);

% set start time as well
params.tStart = 0.5;
params.tFinal = 1;
[t,x] = simulate(sys,params);

% check if start and end time correct
res(end+1,1) = withinTol(t(1),params.tStart) && withinTol(t(end),params.tFinal);

% set time step size as well
params.timeStep = 0.1;
params.tStart = 0.5;
params.tFinal = 1;
steps = round((params.tFinal-params.tStart)/params.timeStep);
[t,x,~,y] = simulate(sys,params);

% check if start and end time correct
res(end+1,1) = withinTol(t(1),params.tStart) ...
    && length(t) == steps + 1 && size(x,1) == steps + 1 && size(y,1) == steps + 1 ...
    && withinTol(t(end),params.tFinal);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
