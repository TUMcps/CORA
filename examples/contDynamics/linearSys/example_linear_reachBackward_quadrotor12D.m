function res = example_linear_reachBackward_quadrotor12D
% example_linear_reachBackward_quadrotor12D - example for backward 
%    reachability analysis
%
% Syntax:
%    res = example_linear_reachBackward_quadrotor12D
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] F. Gruber and M. Althoff, "Scalable Robust Safety Filter with
%        Unknown Disturbance Bounds", IEEE TAC, 2022.
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

% init quadrotor system
[A,B,E] = quadrotor();
sys = linearSys('sys',A,B,[],[],[],[],E);

% target set: obstacle
lb = [2;   2;   2;   0;  0.5; 0.5; -pi/4; -pi/4; -pi/4; -0.1; -0.1; -0.1];
ub = [2.5; 2.5; 2.5; 2;  1;   1;    pi/4;  pi/4;  pi/4;  0.1;  0.1;  0.1];
Xend_interval = interval(lb,ub);
params.R0 = polytope(Xend_interval);

% set of controllable inputs [1, Sec. V-D]
U_interval = interval([-9.81;-0.5;-0.5;-0.5],[2.38;0.5;0.5;0.5]);
params.U = zonotope(U_interval);

% set of uncontrollable disturbances [1, Sec. V-D]
params.W = zonotope(zeros(3,1),0*eye(3));

% time horizon
params.tStart = 0;
params.tFinal = 0.5;

% time step size and type of reachability
options.timeStep = 0.01;
options.linAlg = 'outer:AE:timeinterval';

% backward analysis
R = reachBackward(sys,params,options);

% visualization
% figure; hold on; box on;
% projDims = {[1,2],[1,3],[1,2]};
% for p=1:length(projDims)
%     subplot(2,2,p); hold on; box on;
%     useCORAcolors("CORA:contDynamics");
%     plot(R,projDims{p});
%     plot(Xend_interval,projDims{p},'EdgeColor','k','FaceColor','w');
% end


% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
