function res = test_linearSys_particularSolution_constant
% test_linearSys_particularSolution_constant - unit test for the
%    computation of the particular solution for constant inputs
%
% Syntax:
%    res = test_linearSys_particularSolution_constant
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

% init system, input, and algorithm parameters
A = [-1 -4; 4 -1];
sys = linearSys(A);
u = [1;0];
U = zonotope([1;0],[0.5 1; -1 1]);
timeStep = 0.05;
truncationOrder = 6;

% function for analytical solution (2D): A^-1 (e^Adt - I) * u
true_sol = @(u,t) inv(A)*(expm(A*t) - eye(2)) * u;

% compute particular solution (vector, set)
[Ptp_u,C_input_u,Pti_u] = ...
    particularSolution_constant(sys,u,timeStep,truncationOrder);
[Ptp_U,C_input_U,Pti_U] = ...
    particularSolution_constant(sys,U,timeStep,truncationOrder);

% compare to analytical solution
Ptp_u_true = true_sol(u,timeStep);
Ptp_U_true = true_sol(U,timeStep);
assert(compareMatrices(Ptp_u,Ptp_u_true,tol,"equal",true));
assert(isequal(Ptp_U,Ptp_U_true,tol));

% since u in U, containment for particular solutions follows
assert(contains(Ptp_U,Ptp_u));
assert(contains(Pti_U,Pti_u));

% time-interval solution must contain all true solutions with t in [0,dt],
% u in U; we choose the vertices in U for this
V = vertices(U);
t = linspace(0,timeStep,101);
eAt_V = cell2mat(arrayfun(@(t_) true_sol(V,t_),t,'UniformOutput',false));
assert(all(contains(Pti_U,eAt_V,'exact',tol)));

% the error term is computed with the following meaning: the true input
% solution A^-1(e^At - I)*u with t in [0,dt], u in U must be contained in 
% the linear interpolation between 0 and A^-1(e^At - I)*u at t=dt
% enlarged by the error term
for i=1:100
    % sample state and time
    t = rand*timeStep;
    u = randPoint(U);

    % true propagated input solution
    u_prop_true = true_sol(u,t);

    % linear interpolation
    u_prop_dt = true_sol(u,timeStep);
    u_prop_lint = t/timeStep * u_prop_dt;

    % check containment
    assertLoop(contains(u_prop_lint + C_input_U,u_prop_true,'exact',tol), i);
end


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
