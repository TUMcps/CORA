function res = test_nonlinearARX_computeGO_02_linsys
% test_nonlinearARX_computeGO_02_linsys - unit test for computing a GO
%   model for a nonlinearARX object
%
% Syntax:
%    res = test_nonlinearARX_computeGO_02_linsys
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Laura Luetzow
% Written:       05-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n_k = 10; % number of time steps

% set random number stream
rng('default')

% System Dynamics ---------------------------------------------------------

dt = 0.1; % sampling time
n_u = 1;
n_y = 2;
n_p = 2;

A_bar = {[-0.4 0.6; 0.6 -0.4];[0.1 0; 0.2 -0.5]};
B_bar = {[0; 0];[0.3; -0.7];[0.1; 0]};
sys_LTI = linearARX('exampleSys',A_bar, B_bar,dt);

fun = @(y,u) A_bar{1} * y(3:4) + A_bar{2} * y(1:2) + B_bar{1} * u(3) +...
    B_bar{2} * u(2) + B_bar{3} * u(1);
sys_NL = nonlinearARX('exampleSysNL', fun,dt,n_y, n_u, n_p);

% compute symbolic derivatives
options.tensorOrder = 2;
options.tensorOrderOutput = 2;
derivatives(sys_NL,options);

% GO Models ---------------------------------------------------------------

x0_ref = randn(n_y*n_p,1);
u_ref = randn(n_u,n_k);

p_GO_LTI = computeGO(sys_LTI, x0_ref, u_ref, n_k);
p_GO_NL = computeGO(sys_NL, x0_ref, u_ref, n_k);


% check equality
for k=1:n_k
    assertLoop(sum(abs(p_GO_LTI.A{k} - p_GO_NL.A{k}), 'all') <= 1e-6,k)
    assertLoop(sum(abs(p_GO_LTI.C{k} - p_GO_NL.C{k}), 'all') <= 1e-6,k)

    for j=1:k
        assertLoop(sum(abs(p_GO_LTI.B{k,j} - p_GO_NL.B{k,j}), 'all') <= 1e-6,k,j)
        assertLoop(sum(abs(p_GO_LTI.D{k,j} - p_GO_NL.D{k,j}), 'all') <= 1e-6,k,j)
    end
end

res = true;
end

% ------------------------------ END OF CODE ------------------------------
