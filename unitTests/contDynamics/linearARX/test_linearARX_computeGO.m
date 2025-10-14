function res = test_linearARX_computeGO
% test_linearARX_computeGO - unit test for computing a GO model for a
%   linearARX object
%
% Syntax:
%    res = test_linearARX_computeGO
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
sys = linearARX('exampleSys',A_bar, B_bar,dt);

% GO Models ---------------------------------------------------------------

x0_ref = randn(n_y*n_p,1);
u_ref = randn(n_u,n_k);

p_GO = computeGO(sys, x0_ref, u_ref, n_k);

x0 = randn(n_y*n_p,1);
u = randn(n_u,n_k);
params.x0 = x0;
params.u = u;
params.tFinal = (n_k-1)*dt;

[~,~,~,y_sim] = simulate(sys,params);

% check equality
for k=n_p:n_k-1
    x = p_GO.x(:,k+1) + p_GO.A{k+1} * (x0 - x0_ref);
    y = p_GO.y(:,k+1) + p_GO.C{k+1} * (x0 - x0_ref);
    for i=0:k
        x = x + p_GO.B{k+1,i+1} * (u(:,i+1) - u_ref(:,i+1));
        y = y + p_GO.D{k+1,i+1} * (u(:,i+1) - u_ref(:,i+1));
    end
    
    assert(~(sum(y-y_sim(:,k+1), 'all') > 1e-6 || ...
            sum(x(end-n_p+1:end)-y_sim(:,k+1), 'all') > 1e-6));
end
res = true;
end

% ------------------------------ END OF CODE ------------------------------
