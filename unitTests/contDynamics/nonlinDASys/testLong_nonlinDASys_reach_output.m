function res = testLong_nonlinDASys_reach_output
% testLong_nonlinDASys_reach_output - tests if output equation works
%
% Syntax:
%    res = testLong_nonlinDASys_reach_output
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
params.y0guess = 0;

% reachability settings
options.timeStep = 0.01;
options.taylorTerms = 4;
options.errorOrder = 5;
options.zonotopeOrder = 30;
options.alg = 'lin';
options.tensorOrder = 2;
options.tensorOrderOutput = 2;

% dynamic equation
f = @(x,y,u) [x(1)^2 - x(2) + y(1); x(2) - u(1)];
% constraint equation
g = @(x,y,u) y(1);
% linear output equation
h_lin = @(x,y,u) x(1) + x(2);
% quadratic output equation
h_quad = @(x,y,u) x(1) + x(2)^2;
% cubic output equation
h_cub = @(x,y,u) x(1) + x(2)^3;

% instantiate nonlinearSys objects
sys1 = nonlinDASys('sys1',f,g,h_lin);
sys2 = nonlinDASys('sys2',f,g,h_quad);
sys3 = nonlinDASys('sys3',f,g,h_cub);

% reachability analysis
reach(sys1,params,options);
reach(sys2,params,options);
options.tensorOrderOutput = 3;
reach(sys3,params,options);

% ------------------------------ END OF CODE ------------------------------
