function res = test_linearSys_particularSolution_timeVarying
% test_linearSys_particularSolution_timeVarying - unit test for the
%    computation of the particular solution due to time-varying inputs
%
% Syntax:
%    res = test_linearSys_particularSolution_timeVarying
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
U = zonotope([1;0],[0.5 1; -1 1]);
timeStep = 0.05;
truncationOrder = 6;

% constant input solution in 2D
const_sol = @(u,t) inv(A)*(expm(A*t) - eye(2)) * u;

% compute particular solution
Ptp = particularSolution_timeVarying(sys,U,timeStep,truncationOrder);

% since U contains the origin...
% ...Ptp must contain the origin
assert(contains(Ptp,zeros(2,1),'exact',tol));

% time-point particular solution for time-varying inputs must contain
% time-point particular solution for time-constant inputs
Ptp_const = const_sol(U,timeStep);
assert(contains(Ptp,Ptp_const,'exact',tol));

% time-interval particular solution for time-varying inputs must contain
% particular solution for time-constant inputs at any point in time within
% the given time interval
% for i=1:10
%     t = rand*timeStep;
%     assertLoop(contains(Pti,const_sol(U,t),'exact',tol),i);
% end


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
