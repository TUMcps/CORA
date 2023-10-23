function res = test_linearARMAX_simulate
% test_linearARMAX_simulate - unit test for simulate
%
% Syntax:
%    res = test_linearARMAX_simulate
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
tol_ARMAX = 1e-5; %1e-6 is fine as well
tol_time = 1e-10;

% time horizon
tFinal = 2; % time horizon
dt = 0.1; % or 0.05 or 0.1 % sampling time

res = false;

%% define linearSys -------------------------------------------------------

% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
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

% true initial state
x0 = 10*ones(dim_x,1);

% initialize linearSys object
sys_lin = linearSys(A,B,[],C,D);
sys_linDT = linearSysDT(sys_lin, dt);
sys_ARMAX = linearARMAX(sys_linDT);

% model parameters --------------------------------------------------------

params.tFinal = tFinal;
dt_steps = tFinal / dt;
% initial state
params.x0 = x0; 
% disturbance sets
W = zonotope(0.02+zeros(dim_x,1),0.02*diag(ones(dim_x,1)));
V = zonotope(-0.01+zeros(dim_y,1),0.01*diag(ones(dim_y,1)));


%% simulate ---------------------------------------------------------------

% vectors for u, w, and v
params.u(1,:) = 80*sin(8*[0:dt:tFinal]); 
params.u(2,:) = 20*randn(dim_u-2,dt_steps+1);
params.u(3,:) = 10*sin(5*[0:dt:tFinal]+1); 
params.w  = randPoint(W,dt_steps);
params.v = randPoint(V,dt_steps+1);

% simulate the linear discrete-time system
[t_linDT,~,~,y_linDT] = simulate(sys_linDT,params);

% simulate the linear ARMAX system with normal parametrization
params.u = [params.u; params.w zeros(dim_x,1); params.v];
start_idx = sys_ARMAX.dim +1;
params.y0 = y_linDT(1:start_idx-1, : )';
[t_ARMAX,~,~,y_ARMAX] = simulate(sys_ARMAX,params);

% simulate the linear ARMAX system with time-varying parametrization
setTVP(sys_ARMAX);
[t_tvpARMAX,~,~,y_tvpARMAX] = simulate(sys_ARMAX,params);

%% results
% check if ARMAX and SS give the same simulation results
res_DT_ARMAX = withinTol(y_linDT,y_ARMAX,tol_ARMAX*dt_steps);
res_DT_ARMAX = res_DT_ARMAX & withinTol(t_linDT,t_ARMAX,tol_time);

% check if normal ARMAX formulation and ARMAX with time-varying parameters 
% give the same simulation results
res_ARMAX_tvpARMAX = withinTol(y_ARMAX,y_tvpARMAX,tol_ARMAX*dt_steps);
res_ARMAX_tvpARMAX = res_ARMAX_tvpARMAX & withinTol(t_ARMAX,t_tvpARMAX,tol_time);

% all checks ok
res = all(res_DT_ARMAX & res_ARMAX_tvpARMAX, 'all');

end

% ------------------------------ END OF CODE ------------------------------
