function res = test_linearSys_reachBackward_2D
% test_linearSys_reachBackward_2D - unit test for backward reachability
%    analysis examining the EA and AE case for different input and
%    disturbance capacities
%
% Syntax:
%    res = test_linearSys_reachBackward_2D
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] F. Gruber and M. Althoff, "Computing Safe Sets of Linear
%        Sampled-Data Systems", 2020.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       20-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init 2D quadrotor system
A = [0 1; 0 0];
B = [0; 1];
E = [1 0; 0 1];
sys = linearSys(A,B,[],[],[],[],E);

% model parameters
params.tStart = 0;
params.tFinal = 0.1;
% take target set from [1]
Z = zonotope(zeros(2,1),[-0.3022, 0.0540; 0.1374, 0.1187]);
params.R0 = polytope(Z);

% only time step size required...
options.timeStep = 0.01;

% EA reachability
options.linAlg = 'inner:EA:timepoint';

% backward reachability analysis
% - compute without inputs, with disturbances
params.U = zonotope(0);
params.W = zonotope(zeros(2,1),0.01*eye(2));
R_small = reachBackward(sys,params,options);
% - compute without inputs or disturbances
params.U = zonotope(0);
params.W = zonotope(zeros(2,1));
R_medium = reachBackward(sys,params,options);
% - compute with inputs, without disturbances
params.U = zonotope(0,1);
params.W = zonotope(zeros(2,1));
R_large = reachBackward(sys,params,options);

% the more input capacity vs. disturbance capacity, the larger the backward
% reachable set must be
assert(contains(R_large.timePoint.set{end},R_medium.timePoint.set{end}))
assert(contains(R_medium.timePoint.set{end},R_small.timePoint.set{end}))


% AE reachability
options.linAlg = 'outer:AE:timepoint';

% backward reachability analysis
% - compute with inputs, without disturbances
params.U = zonotope(0,1);
params.W = zonotope(zeros(2,1));
R_small = reachBackward(sys,params,options);
% - compute without inputs or disturbances
params.U = zonotope(0);
params.W = zonotope(zeros(2,1));
R_medium = reachBackward(sys,params,options);
% - compute without inputs, with disturbances
params.U = zonotope(0);
params.W = zonotope(zeros(2,1),0.1*eye(2));
R_large = reachBackward(sys,params,options);

% the more disturbance capacity vs. input capacity, the larger the backward
% reachable set must be
assert(contains(R_large.timePoint.set{end},R_medium.timePoint.set{end}))
assert(contains(R_medium.timePoint.set{end},R_small.timePoint.set{end}))


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
