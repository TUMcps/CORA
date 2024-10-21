function res = test_linearSys_reach_01
% test_linearSys_reach_01 - unit test function of linear reachability
%    analysis; checks the solution of the linearSys class for a small
%    example without inputs against the analytical solution e^At*x0
%
% Syntax:
%    res = test_linearSys_reach_01
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Hendrik Roehm, Matthias Althoff
% Written:       02-March-2016
% Last update:   03-March-2016 (HR)
%                12-August-2016 (MA)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-10;

timeStep = 0.1;
numberOfTimeSteps = 50;

% dynamics
A = [0.1 1; -1 0.1];
linsys = linearSys('linearSys',A);

% model parameters
params.R0 = zonotope([[1; 1], diag([0.5, 0.5])]);
params.U = zonotope(0);
params.uTrans = 0;
params.tFinal = timeStep * numberOfTimeSteps;

options.timeStep = timeStep;
options.taylorTerms = 10;
options.zonotopeOrder = 20;

% compute reachable set
R = reach(linsys,params,options);

% read out vertices of last set
Vcomputed = vertices(R.timePoint.set{end});

% check solution against analytical solution by comparison of the vertices
%    x(t) = e^(A*t)*x0
Vexact = vertices(expm(A)^(timeStep*numberOfTimeSteps)*params.R0);

% compare
assert(compareMatrices(Vexact,Vcomputed,tol,"equal",false));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
