function res = test_linearSys_homogeneousSolution_blocks
% test_linearSys_homogeneousSolution_blocks - unit test for the
%    computation of the homogeneous time-point and time-interval solution
%    using block decomposition
%
% Syntax:
%    res = test_linearSys_homogeneousSolution_blocks
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
tol = 1e-6;

% init algorithm parameters
timeStep = 0.01;
truncationOrder = 5;

% init sparse system
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
sys = linearSys(A,eye(5));
X = zonotope([40;20;10;-20;30],eye(5));

% compute homogeneous solutions without blocks
[Htp,Hti] = homogeneousSolution(sys,X,timeStep,truncationOrder);
% compute homogeneous solutions with blocks
[Htp_blocks,Hti_blocks] = homogeneousSolution(sys,X,timeStep,truncationOrder);
Htp_blocks = compact(recompose(Htp_blocks),'zeros');
Hti_blocks = compact(recompose(Hti_blocks),'zeros');

% time-point solution must be equal
assert(isequal(Htp,Htp_blocks,tol));
% decomposed time-interval solution must contain non-decomposed solution
assert(contains(Hti_blocks,Hti,'exact:polymax',tol));


% dense system matrix
A = [-1 -4 1 -2 1; 4 -1 -1 2 3; 5 2 -3 1 -2; -2 1 -1 -3 1; 1 -1 5 4 -2];
sys = linearSys(A,eye(5));
X = zonotope([40;20;10;-20;30],eye(5));

% compute homogeneous solutions without blocks
[Htp,Hti] = homogeneousSolution(sys,X,timeStep,truncationOrder);
% compute homogeneous solutions with blocks
[Htp_blocks,Hti_blocks] = homogeneousSolution(sys,X,timeStep,truncationOrder);
Htp_blocks = compact(recompose(Htp_blocks),'zeros');
Hti_blocks = compact(recompose(Hti_blocks),'zeros');

% decomposed solutions must contain non-decomposed solutions
assert(contains(Htp_blocks,Htp,'exact',tol));
assert(contains(polytope(Hti_blocks),Hti,'exact',tol));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
