function res = test_nonlinearSysDT_computeGO_02_linsys
% test_nonlinearSysDT_computeGO_02_linsys - unit test for comparing the GO 
%   model of a linearSYSDT object with the GO model of a nonlinearSysDT 
%   object
%
% Syntax:
%    res = test_nonlinearSysDT_computeGO_02_linsys
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Laura Luetzow
% Written:       04-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n_k = 10; % number of time steps

% set random number stream
rng('default')

% System Dynamics ---------------------------------------------------------

dt = 0.1; % sampling time
n_x = 6;
n_u = 4;
n_y = 2;
sys_LTI = linearSysDT('exampleSys',rand(n_x,n_x), rand(n_x,n_u), [], rand(n_y,n_x), rand(n_y,n_u),dt);

fun = @(x,u) sys_LTI.A * x + sys_LTI.B * u;
out_fun = @(x,u) sys_LTI.C * x + sys_LTI.D * u;
sys_NL = nonlinearSysDT('exampleSysNL', fun, dt, n_x, n_u, out_fun, n_y);

% compute symbolic derivatives
options.tensorOrder = 2;
options.tensorOrderOutput = 2;
derivatives(sys_NL,options);

% GO Models ---------------------------------------------------------------

x0_ref = randn(n_x,1);
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
