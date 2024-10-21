function res = test_nonlinearSys_linearize
% test_nonlinearSys_linearize - unit_test_function of linearizing nonlinear 
%    dynamics: Checks the linearization of the nonlinearSys class
%    for the 6D tank example; It is checked whether the A and B matrix
%    are correct for a particular linearization point
%
% Syntax:
%    res = test_nonlinearSys_linearize
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       30-July-2017
% Last update:   12-September-2017
%                13-October-2024 (MW, update syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% tolerance
tol = 1e-12;

% model parameters
params.U = zonotope(0,0.005);
params.uTrans = 0;

% reachability settings
options.alg = 'lin';
options.tensorOrder = 2;

% system dynamics
tank = nonlinearSys(@tank6Eq);
derivatives(tank);

% linearize system
dim_x = 6;
R0 = zonotope([2; 4; 4; 2; 10; 4], 0.2*eye(dim_x));
[~,linsys,linOptions] = linearize(tank,R0,params,options);

% ground truth
A_true = [-0.023490689645049, 0, 0, 0, 0, -0.010000000000000; ...
           0.023490689645049, -0.016610425942763, 0, 0, 0, 0; ...
           0, 0.016610425942763, -0.016610425942763, 0, 0, 0; ...
           0, 0, 0.016610425942763, -0.023490689645049, 0, 0; ...
           0, 0, 0, 0.023490689645049, -0.010505355776936, 0; ...
           0, 0, 0, 0, 0.010505355776936, -0.016610425942763];
U_true = zonotope(zeros(dim_x,1), [0.005; zeros(dim_x-1,1)]);

% compare with obtained values
assert(compareMatrices(linsys.A,A_true,tol,"equal",true));
assert(isequal(linOptions.U,U_true,tol));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
