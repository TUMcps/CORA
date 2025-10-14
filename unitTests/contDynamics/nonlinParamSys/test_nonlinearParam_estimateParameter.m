function completed = test_nonlinearParam_estimateParameter()
% test_nonlinearParam_estimateParameter - unit test for paramter estimation 
%    of nonlinearParam systems.
%
% Syntax:
%    completed = test_nonlinearParam_estimateParameter()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 

% Authors:       Laura Luetzow
% Written:       25-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system dynamics
p_true = 2.8;
f = @(x,u,p) [x(1)*cos(x(2)); ...
    x(1)^2/p(1) * tan(u(1))];
dim_x = 2;
sys = nonlinParamSys(f);

% simulate ground-truth system
simOpts.x0 = [-1; 1];
simOpts.tFinal = 1;
simOpts.u = 0.3;
simOpts.p = p_true;
[t,x] = simulate(sys,simOpts);
u = repmat(simOpts.u,[1,length(t)]);

p_est = estimateParameter(sys,x,t,u);

% simulate system with estimated parameter
simOpts.p = p_est;
[t_est,x_est] = simulate(sys,simOpts);

% check size
assert(size(x,1) == dim_x)
assert(size(x_est,1) == dim_x)

% check deviation
assert(sum(abs(t-t_est),'all') < 1e-9)
assert(sum(abs(x-x_est),'all') < 1e-3)

completed = true;

% ------------------------------ END OF CODE ------------------------------
