function res = testLong_nonlinearSysDT_reach_output
% testLong_nonlinearSysDT_reach_output - tests if output equation works
%
% Syntax:
%    res = testLong_nonlinearSysDT_reach_output
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

% Authors:       Mark Wetzlinger
% Written:       19-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume satisfaction
res = true;

% model parameters
params.tFinal = 0.1;
params.R0 = zonotope(ones(2,1),0.05*eye(2));
params.U = zonotope(0,0.01);

% reachability settings
options.taylorTerms = 4;
options.zonotopeOrder = 30;
options.alg = 'lin';
options.errorOrder = 5;
options.tensorOrder = 2;
options.tensorOrderOutput = 2;

% dynamic equation
f = @(x,u) [x(1)^2 - x(2); x(2) - u(1)];
dt = 0.01;
% linear output equation
g_lin = @(x,u) x(1) + x(2);
% quadratic output equation
g_quad = @(x,u) x(1) + x(2)^2;
% cubic output equation
g_cub = @(x,u) x(1) + x(2)^3;

% instantiate nonlinearSys objects
sys1 = nonlinearSysDT('sys1',f,dt,g_lin);
sys2 = nonlinearSysDT('sys2',f,dt,g_quad);
sys3 = nonlinearSysDT('sys3',f,dt,g_cub);

% reachability analysis
reach(sys1,params,options);
reach(sys2,params,options);
options.tensorOrderOutput = 3;
reach(sys3,params,options);

% ------------------------------ END OF CODE ------------------------------
