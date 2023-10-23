function res = test_nonlinearSysDT_reach_03_linearSysDT
% test_nonlinearSysDT_reach_03_linearSysDT - unit test for reach
%       (comparison with reach from equivalent linearSysDT)
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
U = zonotope([randn(n_u,1) 0.1*eye(n_u)]);
W = zonotope([randn(n_x,1) 0.01*eye(n_x)]);
V = zonotope([randn(n_y,1) 0.02*eye(n_y)]);


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 50;
options.errorOrder = 1;
options.tensorOrder = 2;


% simulate ----------------------------------------------------------------

params_NL = params; 

% simulate linear system
params.uTrans = u_rand;
params.U = U;
R_lin_U = reach(sys_linDT,params,options);

params.W = W;
params.V = V;
R_lin_UWV = reach(sys_linDT,params,options);

% simulate nonlinear system
params_NL.uTrans = [u_rand; zeros(n_x+n_y, dt_steps+1)];
params_NL.U = cartProd(cartProd(U, zonotope(zeros(n_x,1))), zonotope(zeros(n_y,1)));
R_nonlin_U = reach(sys,params_NL,options);

params_NL.U = cartProd(cartProd(U, W), V);
R_nonlin_UWV = reach(sys,params_NL, options);


if isequal(R_lin_U, R_nonlin_U, 0.1) && isequal(R_lin_UWV, R_nonlin_UWV, 0.1)
    % all checks ok
    res = true;
end

end


% ------------------------------ END OF CODE ------------------------------
