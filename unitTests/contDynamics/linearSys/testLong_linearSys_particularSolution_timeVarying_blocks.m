function res = testLong_linearSys_particularSolution_timeVarying_blocks
% testLong_linearSys_particularSolution_timeVarying_blocks - unit test for
%    the computation of the particular solution due to time-varying inputs
%    using block decomposition
%
% Syntax:
%    res = testLong_linearSys_particularSolution_timeVarying_blocks
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

% init input and algorithm parameters
U = 0.5*zonotope([1; 0; 0; 0.5; -0.5],diag([0.2, 0.5, 0.2, 0.5, 0.5]));
timeStep = 0.05;
truncationOrder = 6;

% sparse system matrix
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
sys = linearSys(A,eye(5));
blocks = [1 2; 3 4; 5 5];

% compute particular solution without blocks
Ptp = particularSolution_timeVarying(sys,U,timeStep,truncationOrder);
% compute particular solution with blocks
Ptp_blocks = particularSolution_timeVarying(sys,U,timeStep,truncationOrder,blocks);
Ptp_blocks = compact(recompose(Ptp_blocks),'zeros');

% time-point solutions must be equal
assert(isequal(Ptp,Ptp_blocks,tol));


% dense system matrix
A = [-1 -4 1 -2 1; 4 -1 -1 2 3; 5 2 -3 1 -2; -2 1 -1 -3 1; 1 -1 5 4 -2];
sys = linearSys(A,eye(5));
blocks = [1 2; 3 5];

% compute particular solution without blocks
Ptp = particularSolution_timeVarying(sys,U,timeStep,truncationOrder);
% compute particular solution with blocks
Ptp_blocks = particularSolution_timeVarying(sys,U,timeStep,truncationOrder,blocks);

% compact and reduce to make containment check faster to compute
Ptp_blocks = compact(recompose(Ptp_blocks),'zeros');
Ptp_blocks = reduce(Ptp_blocks,'girard',5);

% time-point solutions must be equal
assert(contains(polytope(Ptp_blocks),Ptp,'exact',tol));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
