function res = example_linear_reachBackward_2D
% example_linear_reachBackward_2D - example for backward reachability
%    analysis, inspired by [1, Sec. V-A]
%
% Syntax:
%    res = example_linear_reachBackward_2D
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
params.tFinal = 3;
% take target set from [1]
Z = zonotope(zeros(2,1),[-0.3022, 0.0540; 0.1374, 0.1187]);
params.R0 = polytope(Z);
params.U = zonotope(0,1);
params.W = zonotope(zeros(2,1),0.1*eye(2));

% reachability settings
options.timeStep = 0.1;
options.linAlg = 'inner:EA:timeinterval';
options.verbose = true;

% backward reachability analysis
R = reachBackward(sys,params,options);

% visualization
figure; hold on; box on;
useCORAcolors("CORA:contDynamics")
plot(R,[1,2]);
plot(Z,[1,2],'k');


% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
