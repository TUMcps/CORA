function res = test_linearSys_oneStep
% test_linearSys_oneStep - unit test for the computation of the one-step
%    reachable set
%
% Syntax:
%    res = test_linearSys_oneStep
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
U = zonotope([1;0],[0.5 1; -1 1]);
u = [0;0];
timeStep = 0.05;
truncationOrder = 6;

% compute reachable sets of first step
[Rtp,Rti,Htp,Hti,PU,Pu,C_state,C_input] = ...
    oneStep(sys,X,U,u,timeStep,truncationOrder);

% since the input set contains the origin, the time-point/interval
% reachable set must contain the time-point/interval homogeneous solution
assert(contains(Rtp,Htp,'exact',tol));
assert(contains(Rti,Hti,'exact',tol));

% since we have no input vector, the particular solution due to the
% constant offset vector and its error term must be zero
assert(compareMatrices(Pu,zeros(2,1),tol,'equal',true));
assert(all(C_input == zeros(2,1)));

% the affine time-point solution is e^At*X + Pu, so here: e^At*X
assert(isequal(Htp,expm(A*timeStep)*X,tol));

% the affine time-interval solution is enclose(X,Htp) + C_state
assert(isequal(Hti,enclose(X,expm(A*timeStep)*X)+C_state,tol));

% the time-point/interval reachable set is the Minkowski sum of the
% time-point/interval affine solution and the particular solution 
assert(isequal(Rtp,Htp+PU,tol));
assert(isequal(Rti,Hti+PU,tol));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
