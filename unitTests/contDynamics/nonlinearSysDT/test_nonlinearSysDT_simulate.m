function res = test_nonlinearSysDT_simulate
% test_nonlinearSysDT_simulate - unit test for simulate
%
% Syntax:
%    res = test_nonlinearSysDT_simulate
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
% Written:       23-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = false;

% define linearSysDT ------------------------------------------------------

% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
n_x = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];
n_u = size(B,2);

% output matrix: q x n
C = [1 1 0 0;
     0 -0.5 0.5 0];
n_y = 2;

% feedthrough matrix: q x m
D = [0 0 1;
     0 0 0];

% initialize linearSysDT-object
sys_lin = linearSys(A,B,[],C,D);
dt = 0.05;
sys_linDT = linearSysDT(sys_lin,dt);

% initialize nonlinearSysDT-object
B_new = [sys_linDT.B eye(n_x) zeros(n_x,n_y)]; % from combining inputs and disturbances in u
D_new = [sys_linDT.D zeros(n_y,n_x), eye(n_y)]; % from combining inputs and disturbances in u
fun = @(x,u) sys_linDT.A*x + B_new * u;
out_fun = @(x,u) sys_linDT.C*x + D_new * u;
sys = nonlinearSysDT(fun, dt, n_x, n_u+n_x+n_y, out_fun, n_y);

% model parameters --------------------------------------------------------

% time horizon
params.tFinal = 1;
dt_steps = params.tFinal / dt;
% initial set
params.R0 = zonotope(10*ones(n_x,1),0.5*diag(ones(n_x,1)));
% initial state
params.x0 = center(params.R0);
% vectors for u, w, and v
u_rand = randn(n_u, dt_steps+1);
w_rand = randn(n_x, dt_steps);
v_rand = randn(n_y, dt_steps+1);


% simulate ----------------------------------------------------------------

params_NL = params; 

% simulate linear system
params.u = u_rand;
[t_lin_u, x_lin_u, ~, y_lin_u] = simulate(sys_linDT,params);

params.w = w_rand;
params.v = v_rand;
[t_lin_uwv, x_lin_uwv, ~, y_lin_uwv] = simulate(sys_linDT,params);

% simulate nonlinear system
params_NL.u = [u_rand; zeros(n_x+n_y, dt_steps+1)];
[t_nonlin_u, x_nonlin_u, ~, y_nonlin_u] = simulate(sys,params_NL);

params_NL.u = [u_rand; [w_rand zeros(n_x, 1)]; v_rand];
[t_nonlin_uwv, x_nonlin_uwv, ~, y_nonlin_uwv] = simulate(sys,params_NL);


if isequal(t_lin_u, t_nonlin_u, t_lin_uwv, t_nonlin_uwv) ...
        && sum(abs(y_lin_u - y_nonlin_u),'all') < 1e-6 ...
        && sum(abs(y_lin_uwv - y_nonlin_uwv),'all') < 1e-6 ...
        && sum(abs(x_lin_u - x_nonlin_u),'all') < 1e-6 ...
        && sum(abs(x_lin_uwv - x_nonlin_uwv),'all') < 1e-6
    % all checks ok
    res = true;
end

end


% ------------------------------ END OF CODE ------------------------------
