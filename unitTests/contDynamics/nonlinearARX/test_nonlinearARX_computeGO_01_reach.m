function res = test_nonlinearARX_computeGO_01_reach
% test_nonlinearARX_computeGO_01_reach - unit test for computing a GO model for a
%   nonlinearARX object
%
% Syntax:
%    res = test_nonlinearARX_computeGO_01_reach
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Laura Luetzow
% Written:       14-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n_k = 10; % number of time steps

% set random number stream
rng('default')

% System Dynamics ---------------------------------------------------------

f = @(y,u) [0.5*y(1,1) + u(1,1) - cos(u(2,1)); ...
    0.4*y(3,1) + u(2,1)*cos(y(1,1)); 0.6*y(5,1) + u(4,1)*sin(y(1,1))];
dt = 0.1;
dim_y = 3;
dim_u = 2;
p = 2;
sys = nonlinearARX(f,dt,dim_y, dim_u, p);


% Parameters --------------------------------------------------------------

% time horizon

params.tStart = 0;
params.tFinal = dt * (n_k-1);

% initilization
x0_ref = 4*randn(dim_y*p,1);
x0 = randn(dim_y*p,1);

% input
u = 0.1*rand(2,n_k);
u_ref = rand(2,n_k);

% options
options.zonotopeOrder = 50; %zonotope order
options.tensorOrder = 2;
options.errorOrder = 1;
options.lagrangeRem.simplify = 'simplify';


% Linearization -----------------------------------------------------------

% simulate y_nom
params.x0 = x0_ref;
params.u = u_ref;
[~,~,~,y_sim] = simulate(sys, params);
y_sim = y_sim';

% stack the input sets
uref_stacked = getStackedU(sys, u_ref);
u_stacked = getStackedU(sys, u);
derivatives(sys,options);

% compute linearized parameters and nominal solution
p_GO = computeGO(sys, x0_ref, u_ref, n_k);
y_nom = p_GO.y;

% compare y_nom und y_sim
assert(all(y_nom-y_sim < 1e-9, 'all'))


% Simluation without considering Linearization Error ----------------------

x_lin = zeros(p*dim_y, n_k);
y_lin = zeros(dim_y, n_k);
y_lin(:,1:p) = reshape(x0,dim_y,[]);
y_rec = y_lin;
y_lin2 = y_lin;

C_tilde = [zeros(dim_y,(p-1)*dim_y) eye(dim_y)];

for k=p+1:n_k
    xd_last_rec = reshape(y_rec(:,k-p:k-1),[],1);
    xnom_last = reshape(y_nom(:,k-p:k-1),[],1);

    % output using the Jacobians and the previous output
    [A_jac,B_jac] = sys.jacobian(xnom_last, uref_stacked(:,k));
    y_rec(:,k) = A_jac * xd_last_rec + B_jac * u_stacked(:,k);

    % output using the linearized system matrices and the initial outputs
    x_lin(:,k) = p_GO.A{k} * x0;
    for j = 1:k
        x_lin(:,k) = x_lin(:,k) + p_GO.B{k,j} * u(:,j);
    end
    y_lin(:,k) = C_tilde * x_lin(:,k);

    % output using the linearized system matrices and the initial outputs
    y_lin2(:,k) = p_GO.C{k} * x0;
    for j = 1:k
        y_lin2(:,k) = y_lin2(:,k) + p_GO.D{k,j} * u(:,j);
    end
end

% compare outputs
assert(all(y_lin-y_rec < 1e-9, 'all'))
assert(all(y_lin-y_lin2 < 1e-9, 'all'))


% Simluation considering random linearization error -----------------------

x_lin = zeros(p*dim_y, n_k);
y_lin = zeros(dim_y, n_k);
y_lin(:,1:p) = reshape(x0,dim_y,[]);
y_rec = y_lin;
y_lin2 = y_lin;
l = 0.1 * randn(dim_y, n_k);

C_tilde = [zeros(dim_y,(p-1)*dim_y) eye(dim_y)];

for k=p+1:n_k
    xd_last_rec = reshape(y_rec(:,k-p:k-1),[],1);
    xnom_last = reshape(y_nom(:,k-p:k-1),[],1);

    % output using the Jacobians and the previous output
    [A_jac,B_jac] = sys.jacobian(xnom_last, uref_stacked(:,k));
    y_rec(:,k) = A_jac * xd_last_rec + B_jac * u_stacked(:,k) + l(:,k);

    % output using the linearized system matrices and the initial outputs
    x_lin(:,k) = p_GO.A{k} * x0;
    for j = 1:k
        x_lin(:,k) = x_lin(:,k) + p_GO.B{k,j} * u(:,j);
        if j >= p+1
            x_lin(:,k) = x_lin(:,k) + p_GO.F{k,j} * l(:,j);
        end
    end
    y_lin(:,k) = C_tilde * x_lin(:,k);

    % output using the linearized system matrices and the initial outputs
    y_lin2(:,k) = p_GO.C{k} * x0;
    for j = 1:k
        y_lin2(:,k) = y_lin2(:,k) + p_GO.D{k,j} * u(:,j);
        if j >= p+1
            y_lin2(:,k) = y_lin2(:,k) + p_GO.E{k,j} * l(:,j);
        end
    end
end

% compare yd_lin und yd_rec
assert(all(y_lin-y_rec < 1e-9, 'all'))
assert(all(y_lin-y_lin2 < 1e-9, 'all'))


res = true;
end

% ------------------------------ END OF CODE ------------------------------
