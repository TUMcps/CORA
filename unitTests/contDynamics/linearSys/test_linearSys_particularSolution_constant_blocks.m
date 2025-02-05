function res = test_linearSys_particularSolution_constant_blocks
% test_linearSys_particularSolution_constant_blocks - unit test for the
%    computation of the particular solution for constant inputs
%
% Syntax:
%    res = test_linearSys_particularSolution_constant_blocks
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
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-7;

% init system, input, and algorithm parameters
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
sys = linearSys(A,eye(5));
U_vector = [1; 2; 0; -3; 1];
U_set = 0.5*zonotope([1; 0; 0; 0.5; -0.5],diag([0.2, 0.5, 0.2, 0.5, 0.5]));
timeStep = 0.05;
truncationOrder = 6;
blocks = [1 2; 3 4; 5 5];

% function for analytical solution (5D): A^-1 (e^Adt - I) * u
true_sol = @(A,u,t) inv(A)*(expm(A*t) - eye(5)) * u;

% evaluate analytical solution
Pu_true = true_sol(A,U_set,timeStep);
% compute particular solution without blocks
Ptp_u = particularSolution_constant(sys,U_set,timeStep,truncationOrder);
% compute particular solution with blocks and simplify
Ptp_u_blocks = particularSolution_constant(sys,U_set,timeStep,truncationOrder,blocks);
Ptp_u_blocks = compact(recompose(Ptp_u_blocks),'zeros');

% ensure that all are equal (blocks are independent from one another)
assert(isequal(Pu_true,Ptp_u,tol));
assert(isequal(Pu_true,Ptp_u_blocks,tol));

% ...same for vector
Pu_true = true_sol(A,U_vector,timeStep);
Ptp_u = particularSolution_constant(sys,U_vector,timeStep,truncationOrder);
Ptp_u_blocks = particularSolution_constant(sys,U_vector,timeStep,truncationOrder,blocks);
Ptp_u_blocks = recompose(Ptp_u_blocks);
assert(compareMatrices(Pu_true,Ptp_u,tol,"equal",true));
assert(compareMatrices(Pu_true,Ptp_u_blocks,tol,"equal",true));


% dense system matrix (still invertible)
A = [-1 -4 1 -2 1; 4 -1 -1 2 3; 5 2 -3 1 -2; -2 1 -1 -3 1; 1 -1 5 4 -2];
sys = linearSys(A,eye(5));
U_vector = [1; 2; 0; -3; 1];
U_set = 0.5*zonotope([1; 0; 0; 0.5; -0.5],diag([0.2, 0.5, 0.2, 0.5, 0.5]));
timeStep = 0.05;
truncationOrder = 6;
blocks = [1 2; 3 5];

% compute particular solution
Ptp_u = particularSolution_constant(sys,U_set,timeStep,truncationOrder);
Ptp_u_blocks = particularSolution_constant(sys,U_set,timeStep,truncationOrder,blocks);
Ptp_u_blocks = compact(recompose(Ptp_u_blocks),'zeros');
% decomposed solution must contain non-decomposed solution
assert(contains(Ptp_u_blocks,Ptp_u,'exact',tol));

% compute particular solutions
Ptp_u = particularSolution_constant(sys,U_vector,timeStep,truncationOrder);
Ptp_u_blocks = particularSolution_constant(sys,U_vector,timeStep,truncationOrder,blocks);
Ptp_u_blocks = recompose(Ptp_u_blocks);
% decomposed solution is equal to non-decomposed solution in vector case
assert(compareMatrices(Ptp_u,Ptp_u_blocks,tol,"equal",true));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
