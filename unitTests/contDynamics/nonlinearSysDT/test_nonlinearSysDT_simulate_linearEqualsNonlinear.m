function res = test_nonlinearSysDT_simulate_linearEqualsNonlinear
% test_nonlinearSysDT_simulate_linearEqualsNonlinear - unit test to check
%    whether the simulation of a nonlinear discrete-time system with linear
%    dynamics gives the same result as a linear discrete-time system
%
% Syntax:
%    res = test_nonlinearSysDT_simulate_linearEqualsNonlinear
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       22-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
tol = 1e-12;

% time step size
dt = 0.1;

% init linear system
A = [0.84 -0.01; -0.02 0.85];
B = [0.06; -0.01];
C = [0.5 -0.25];
D = 2;
sys_lin = linearSysDT(A,B,dt);
sys_lin_withoutput = linearSysDT(A,B,[],C,D,[],dt);

% init equivalent nonlinear system
f = @(x,u) [0.84*x(1) - 0.01*x(2) + 0.06*u(1);
            -0.02*x(1) + 0.85*x(2) - 0.01*u(1)];
g = @(x,u) 0.5*x(1) - 0.25*x(2) + 2*u(1);
sys_nonlin = nonlinearSysDT(f,dt);
sys_nonlin_withoutput = nonlinearSysDT(f,dt,g);


% model parameters (input varies, see below)
params.tStart = 0;
params.tFinal = 1;
params.x0 = [5;-3];

% no output equation
params.u = 0;
[~,x_lin] = simulate(sys_lin,params);
[~,x_nonlin] = simulate(sys_nonlin,params);
assert(compareMatrices(x_lin,x_nonlin,tol,'equal',true));

% with output equation
params.u = 0;
[~,~,~,y_lin] = simulate(sys_lin_withoutput,params);
[~,~,~,y_nonlin] = simulate(sys_nonlin_withoutput,params);
assert(compareMatrices(y_lin,y_nonlin,tol,'equal',true));

% no output equation
params.u = 2;
[~,x_lin] = simulate(sys_lin,params);
[~,x_nonlin] = simulate(sys_nonlin,params);
assert(compareMatrices(x_lin,x_nonlin,tol,'equal',true));

% with output equation
params.u = 2;
[~,~,~,y_lin] = simulate(sys_lin_withoutput,params);
[~,~,~,y_nonlin] = simulate(sys_nonlin_withoutput,params);
assert(compareMatrices(y_lin,y_nonlin,tol,'equal',true));

% no output equation
params.u = [1 1.5 2 0.5 0.25 0 -0.5 -1 -1.5 -0.5];
[~,x_lin] = simulate(sys_lin,params);
[~,x_nonlin] = simulate(sys_nonlin,params);
assert(compareMatrices(x_lin,x_nonlin,tol,'equal',true));

% with output equation
params.u = [1 1.5 2 0.5 0.25 0 -0.5 -1 -1.5 -0.5 0.5];
[~,~,~,y_lin] = simulate(sys_lin_withoutput,params);
[~,~,~,y_nonlin] = simulate(sys_nonlin_withoutput,params);
assert(compareMatrices(y_lin,y_nonlin,tol,'equal',true));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
