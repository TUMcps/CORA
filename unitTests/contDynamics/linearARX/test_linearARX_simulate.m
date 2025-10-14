function res = test_linearARX_simulate
% test_linearARX_simulate - unit test for simulate
%
% Syntax:
%    res = test_linearARX_simulate
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Example:
%
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in
%       Proc. of the 62th International Conference on Decision and Control,
%       2023.

% Authors:       Laura Luetzow
% Written:       06-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set random number stream
rng('default')

% toleranz for equality comparison
tol_ARX = 1e-5; %1e-6 is fine as well
tol_time = 1e-10;

% time horizon
tFinal = 2; % time horizon
dt = 0.1; % or 0.05 or 0.1 % sampling time

% stable system matrix: n x n
A = [0.9618    0.0286    0.0509   -0.0265
    0.0179    1.0257    0.0548    0.0836
    0.0066   -0.0111    0.9378    0.0365
    0.0990   -0.0478    0.0811    0.9747];
dim_x = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
    2 1 0;
    0 0 1;
    0 -2 1];
dim_u = size(B,2);

% output matrix: q x n
C = [1 1 0 0;
    0 -0.5 0.5 0];

% throughput matrix: q x m
D = [0 0 1;
    0 0 0];
dim_y = size(D,1);

% disturbance matrix: n x r
E = eye(dim_x);

% noise matrix: q x s
F = eye(dim_y);

% initialize linearSys object
sys_linDT = linearSysDT(A,B,[],C,D,[],E,F,dt);
sys_ARX = linearARX(sys_linDT);

% model parameters
params.tFinal = tFinal;
dt_steps = tFinal / dt;
params.x0 = 10*ones(dim_x,1); 
W = zonotope(0.02+zeros(dim_x,1),0.02*diag(ones(dim_x,1)));
V = zonotope(-0.01+zeros(dim_y,1),0.01*diag(ones(dim_y,1)));


% simulate different systems and compare results
% vectors for u, w, and v
params.u(1,:) = 80*sin(8*(0:dt:tFinal)); 
params.u(2,:) = 20*randn(dim_u-2,dt_steps+1);
params.u(3,:) = 10*sin(5*(0:dt:tFinal)+1); 
params.w  = randPoint(W,dt_steps);
params.v = randPoint(V,dt_steps+1);

% simulate the linear discrete-time system
[t_linDT,~,~,y_linDT] = simulate(sys_linDT,params);

% simulate the linear ARX system with normal parametrization
params.u = [params.u; params.w zeros(dim_x,1); params.v];
params.y0 = y_linDT(:,1:sys_ARX.n_p);
[t_ARX,~,~,y_ARX] = simulate(sys_ARX,params);

% simulate the linear ARX system with time-varying parametrization
setTVP(sys_ARX);
[t_tvpARX,~,~,y_tvpARX] = simulate(sys_ARX,params);


% check if ARX and SS give the same simulation results
assert(compareMatrices(y_linDT,y_ARX,tol_ARX*dt_steps,'equal',true));
assert(all(withinTol(t_linDT,t_ARX,tol_time)));

% check if normal ARX formulation and ARX with time-varying parameters 
% give the same simulation results
assert(compareMatrices(y_ARX',y_tvpARX',tol_ARX*dt_steps,'equal',true));
assert(all(withinTol(t_ARX',t_tvpARX',tol_time)));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
