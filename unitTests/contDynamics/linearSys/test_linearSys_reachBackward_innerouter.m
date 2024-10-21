function res = test_linearSys_reachBackward_innerouter
% test_linearSys_reachBackward_innerouter - unit test for backward
%    reachability analysis
%
% Syntax:
%    res = test_linearSys_reachBackward_innerouter
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

% tolerance
tol = 1e-10;

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

% only need time step size
options.timeStep = 0.01;

% AE reachability
% - compute without inputs or disturbances
params.U = zonotope(0);
params.W = zonotope(zeros(2,1));
options.linAlg = 'inner:AE:timepoint';
R_inner = reachBackward(sys,params,options);
options.linAlg = 'outer:AE:timepoint';
R_outer = reachBackward(sys,params,options);
assert(contains(R_outer.timePoint.set{end},R_inner.timePoint.set{end},'exact',tol))

% - compute with inputs, without disturbances
params.U = zonotope(0,1);
params.W = zonotope(zeros(2,1));
options.linAlg = 'inner:AE:timepoint';
R_inner = reachBackward(sys,params,options);
options.linAlg = 'outer:AE:timepoint';
R_outer = reachBackward(sys,params,options);
assert(contains(R_outer.timePoint.set{end},R_inner.timePoint.set{end},'exact',tol))

% - compute without inputs, with disturbances
params.U = zonotope(0);
params.W = zonotope(zeros(2,1),0.1*eye(2));
options.linAlg = 'inner:AE:timepoint';
R_inner = reachBackward(sys,params,options);
options.linAlg = 'outer:AE:timepoint';
R_outer = reachBackward(sys,params,options);
assert(contains(R_outer.timePoint.set{end},R_inner.timePoint.set{end},'exact',tol))

% - compute with inputs, with disturbances
params.U = zonotope(0,1);
params.W = zonotope(zeros(2,1),0.1*eye(2));
options.linAlg = 'inner:AE:timepoint';
R_inner = reachBackward(sys,params,options);
options.linAlg = 'outer:AE:timepoint';
R_outer = reachBackward(sys,params,options);
assert(contains(R_outer.timePoint.set{end},R_inner.timePoint.set{end},'exact',tol))


% EA reachability
% - compute without inputs or disturbances
params.U = zonotope(0);
params.W = zonotope(zeros(2,1));
options.linAlg = 'inner:EA:timepoint';
R_inner = reachBackward(sys,params,options);
options.linAlg = 'outer:EA:timepoint';
R_outer = reachBackward(sys,params,options);
assert(contains(R_outer.timePoint.set{end},R_inner.timePoint.set{end},'exact',tol))

% - compute with inputs, without disturbances
params.U = zonotope(0,1);
params.W = zonotope(zeros(2,1));
options.linAlg = 'inner:EA:timepoint';
R_inner = reachBackward(sys,params,options);
options.linAlg = 'outer:EA:timepoint';
R_outer = reachBackward(sys,params,options);
assert(contains(R_outer.timePoint.set{end},R_inner.timePoint.set{end},'exact',tol))

% - compute without inputs, with disturbances
params.U = zonotope(0);
params.W = zonotope(zeros(2,1),0.01*eye(2));
options.linAlg = 'inner:EA:timepoint';
R_inner = reachBackward(sys,params,options);
options.linAlg = 'outer:EA:timepoint';
R_outer = reachBackward(sys,params,options);
assert(contains(R_outer.timePoint.set{end},R_inner.timePoint.set{end},'exact',tol))

% - compute with inputs, with disturbances
params.U = zonotope(0,1);
params.W = zonotope(zeros(2,1),0.01*eye(2));
options.linAlg = 'inner:EA:timepoint';
R_inner = reachBackward(sys,params,options);
options.linAlg = 'outer:EA:timepoint';
R_outer = reachBackward(sys,params,options);
assert(contains(R_outer.timePoint.set{end},R_inner.timePoint.set{end},'exact',tol))


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
