function res = test_nonlinearSysDT_reach_03_linearSysDT
% test_nonlinearSysDT_reach_03_linearSysDT - unit test for reach
%    (comparison with reach from equivalent linearSysDT)
%
% Syntax:
%    res = test_nonlinearSysDT_reach_03_linearSysDT
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
% Written:       26-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% define linearSysDT ------------------------------------------------------

% stable system matrix: n x n (eigenvalues in unit circle)
A = [0.9810    0.0143    0.0262   -0.0140
    0.0079    1.0133    0.0267    0.0416
    0.0029   -0.0054    0.9680    0.0188
    0.0503   -0.0242    0.0411    0.9877];
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

% disturbance matrix: n x r
E = eye(n_x);

% noise matrix: q x s
F = eye(n_y);

% initialize linearSysDT-object
dt = 0.05;
sys_linDT = linearSysDT(A,B,[],C,D,[],E,F,dt);

% initialize nonlinearSysDT-object
B_new = [sys_linDT.B eye(n_x) zeros(n_x,n_y)];
D_new = [sys_linDT.D zeros(n_y,n_x), eye(n_y)];
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
U = zonotope([randn(n_u,1) 0.1*eye(n_u)]);
W = zonotope([randn(n_x,1) 0.01*eye(n_x)]);
V = zonotope([randn(n_y,1) 0.02*eye(n_y)]);


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 50;
options.errorOrder = 1;
options.tensorOrder = 2;


% reachability ------------------------------------------------------------

params_NL = params; 

% reachability for linear system
params.uTrans = u_rand;
params.U = U;
R_lin_U = reach(sys_linDT,params,options);

params.W = W;
params.V = V;
R_lin_UWV = reach(sys_linDT,params,options);

% reachability for nonlinear system
params_NL.uTrans = [u_rand; zeros(n_x+n_y, dt_steps+1)];
params_NL.U = cartProd(cartProd(U, zonotope(zeros(n_x,1))), zonotope(zeros(n_y,1)));
R_nonlin_U = reach(sys,params_NL,options);

params_NL.U = cartProd(cartProd(U, W), V);
R_nonlin_UWV = reach(sys,params_NL,options);

% compare computed reachable sets
assert(isequal(R_lin_U, R_nonlin_U, 0.1));
assert(isequal(R_lin_UWV, R_nonlin_UWV, 0.1));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
