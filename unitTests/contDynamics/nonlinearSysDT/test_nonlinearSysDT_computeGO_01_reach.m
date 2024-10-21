function res = test_nonlinearSysDT_computeGO_01_reach
% test_nonlinearSysDT_computeGO_01_reach - unit test for computing a GO model for a
%   nonlinearSysDT object
%
% Syntax:
%    res = test_nonlinearSysDT_computeGO_01_reach
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Laura Luetzow
% Written:       21-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n_k = 10; % number of time steps

% set random number stream
rng('default')

% System Dynamics ---------------------------------------------------------

dt = 0.015; % sampling time

fun = @(x,u) tank6EqDT(x,u,dt);

n_x = 6;
n_u = 4;
out_fun = @(x,u) [x(1)^2; x(2)] + [sin(u(3)); u(4)];
n_y = 2;
sys = nonlinearSysDT('tank6', fun, dt, n_x, n_u, out_fun, n_y);

% Reachability Analysis ---------------------------------------------------

params.R0 = zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(6)]);
params.U = 0.001*zonotope(randn(n_u,n_u+2));
params.uTransVec = 0.001*randn(n_u,n_k);
params.tFinal = sys.dt * n_k - sys.dt;

options.zonotopeOrder = 200;
options.tensorOrder = 2;
options.errorOrder = 1;
options.tensorOrderOutput = 2;
sys = setHessian(sys,'standard');
sys = setOutHessian(sys,'standard');
[R, ~, Verror] = reach(sys,params,options);
Y_reach = R.timePoint.set;
 
U = cell(n_k,1);
for k=1:n_k
    U{k} = params.U + params.uTransVec(:,k);
end
X0 = params.R0;
u_ref = params.uTransVec;

% Linearization -----------------------------------------------------------

% compute symbolic derivatives
options.tensorOrder = 2;
options.tensorOrderOutput = 2;
derivatives(sys,options);

p_GO = computeGO(sys, center(X0), u_ref, n_k);

% Reachability without considering linearization error --------------------

% initialization
X_rec = cell(n_k+1,1);
Y_rec = cell(n_k,1);
X_lin = cell(n_k+1,1);
Y_lin = cell(n_k,1);

X_rec{1} = X0;
X_lin{1} = X0;
for k=1:n_k
    % reachability using the Jacobians and the previous state set
    [A_jac,B_jac] = sys.jacobian(p_GO.x(:,k), u_ref(:,k));
    [C_jac,D_jac] = sys.out_jacobian(p_GO.x(:,k), u_ref(:,k));
    Y_rec{k} = p_GO.y(:,k) + C_jac * (X_rec{k}-p_GO.x(:,k)) ...
        + D_jac * (U{k}-u_ref(:,k));
    X_rec{k+1} = p_GO.x(:,k+1) + A_jac * (X_rec{k}-p_GO.x(:,k)) ...
        + B_jac * (U{k}-u_ref(:,k));

    % reachability using the linearized system matrices and the initial state set
    Y_lin{k} = p_GO.y(:,k) + p_GO.C{k} * (X0-p_GO.x(:,1));
    X_lin{k+1} = p_GO.x(:,k+1) + p_GO.A{k} * (X0-p_GO.x(:,1));
    for j = 1:k
        Y_lin{k} = Y_lin{k} + p_GO.D{k,j} * (U{j}-u_ref(:,j));
        X_lin{k+1} = X_lin{k+1} + p_GO.B{k,j} * (U{j}-u_ref(:,j));
    end

    assertLoop(isequal(Y_rec{k}, Y_lin{k}, 1e-6),k)
    assertLoop(isequal(X_rec{k+1}, X_lin{k+1}, 1e-6),k)
end

% Reachability considering linearization error ----------------------------

% Initialization
X_rec = cell(n_k+1,1);
Y_rec = cell(n_k,1);
X_lin = cell(n_k+1,1);
Y_lin = cell(n_k,1);

X_rec{1} = X0;
X_lin{1} = X0;

L_xy = cell(n_k,1);

for k=1:n_k
    % reachability using the Jacobians and the previous state set
    [A_jac,B_jac] = sys.jacobian(p_GO.x(:,k), u_ref(:,k));
    [C_jac,D_jac] = sys.out_jacobian(p_GO.x(:,k), u_ref(:,k));
    Y_rec{k} = p_GO.y(:,k) + C_jac * (X_rec{k}-p_GO.x(:,k)) ...
        + D_jac * (U{k}-u_ref(:,k)) + Verror.L_y{k};
    X_rec{k+1} = p_GO.x(:,k+1) + A_jac * (X_rec{k}-p_GO.x(:,k)) ...
        + B_jac * (U{k}-u_ref(:,k)) + Verror.L_x{k};

    % reachability using the linearized system matrices and the initial state set
    Y_lin{k} = p_GO.y(:,k) + p_GO.C{k} * (X0-p_GO.x(:,1));
    X_lin{k+1} = p_GO.x(:,k+1) + p_GO.A{k} * (X0-p_GO.x(:,1));
    L_xy{k} = cartProd(Verror.L_x{k}, Verror.L_y{k});
    for j = 1:k
        Y_lin{k} = Y_lin{k} + p_GO.D{k,j} * (U{j}-u_ref(:,j)) ...
            + p_GO.E{k,j} * L_xy{j};
        X_lin{k+1} = X_lin{k+1} + p_GO.B{k,j} * (U{j}-u_ref(:,j)) ...
            + p_GO.F{k,j} * L_xy{j};
    end

    assertLoop(isequal(Y_rec{k}, Y_lin{k}, 1e-6),k)
    assertLoop(isequal(X_rec{k+1}, X_lin{k+1}, 1e-6),k)

    assertLoop(isequal(compact(Y_reach{k},'aligned',1e-2),compact(Y_lin{k},'aligned',1e-2), 1e-2),k)
end
res = true;
end

% ------------------------------ END OF CODE ------------------------------
