function res = test_linearSysDT_computeGO_01_sim
% test_linearSysDT_computeGO_01_sim - unit test for computing a GO model for a
%   linearARX object
%
% Syntax:
%    res = test_linearSysDT_computeGO_01_sim
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
n_x = 3;
n_u = 1;
n_y = 2;

sys = linearSysDT('exampleSys',rand(n_x,n_x), rand(n_x,n_u), [], rand(n_y,n_x), rand(n_y,n_u),dt);

% GO Models ---------------------------------------------------------------

x0_ref = randn(n_x,1);
u_ref = randn(n_u,n_k);

p_GO = computeGO(sys, x0_ref, u_ref, n_k);

params.x0 = randn(n_x,1);
params.u = randn(n_u,n_k);
params.tFinal = (n_k-1)*dt;

[~,x_sim,~,y_sim] = simulate(sys,params);

% check equality
for k=0:n_k-1
    x = p_GO.x(k+1) + p_GO.A{k+1} * (params.x0 - x0_ref);
    y = p_GO.y(:,k+1) + p_GO.C{k+1} * (params.x0 - x0_ref);
    for i=0:k
        x = x + p_GO.B{k+1,i+1} * (params.u(:,i+1) - u_ref(:,i+1));
        y = y + p_GO.D{k+1,i+1} * (params.u(:,i+1) - u_ref(:,i+1));
    end
    assertLoop(sum(y-y_sim(k+1,:)', 'all') <= 1e-6,i)
    assertLoop(sum(x-x_sim(k+1,:)', 'all') <= 1e-6,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
