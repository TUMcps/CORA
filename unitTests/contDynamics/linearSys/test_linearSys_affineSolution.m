function res = test_linearSys_affineSolution
% test_linearSys_affineSolution - unit test for the computation of the
%    affine solution
%
% Syntax:
%    res = test_linearSys_affineSolution
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
tol = 1e-14;

% init system, state, input, and algorithm parameters
A = [-1 -4; 4 -1];
sys = linearSys(A);
X = zonotope([40;20],[1 4 2; -1 3 5]);
u = [2;-1];
timeStep = 0.05;
truncationOrder = 6;

% compute reachable sets of first step
[Htp,Pu,Hti,C_state,C_input] = ...
    affineSolution(sys,X,u,timeStep,truncationOrder);

% compare particular solution to analytical solution
Pu_true = inv(A)*(expm(A*timeStep) - eye(2)) * u;
assert(compareMatrices(Pu,Pu_true,tol,'equal',true));

% time-interval affine solution must contain time-point affine solution
assert(contains(Hti,Htp,'exact',tol));

% the affine time-point solution is e^At*X + Pu
assert(isequal(Htp,expm(A*timeStep)*X+Pu,tol));

% the affine time-interval solution is enclose(X,Htp+Pu) + error terms
assert(isequal(Hti,enclose(X,expm(A*timeStep)*X+Pu)+C_state+C_input,tol));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
