function res = example_linear_reachBackward_pursuitevasion
% example_linear_reachBackward_pursuitevasion - example for backward
%    reachability analysis using the pursuit evasion game: the target set
%    defines the set of collisions; in the EA case, we compute all states
%    from which player 1 (input) can catch player 2 (disturbance), whereas
%    the AE case computes all states from which player 1 cannot avoid
%    collision
%
% Syntax:
%    res = example_linear_reachBackward_pursuitevasion
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    -
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

% system dynamics: relative position and velocity of two quadrotors
A = [0 1 0 0;
     0 0 0 0;
     0 0 0 1;
     0 0 0 0];
% input (player 1): relative acceleration in x and y
B = [0 0;
     1 0;
     0 0;
     0 1];
% disturbance (player 2): relative acceleration in x and y
E = [0  0;
    -1  0;
     0  0;
     0 -1];
sys = linearSys('sys',A,B,[],[],[],[],E);

% potential collision at time 0
params.tStart = 0;
% backward in time until...
params.tFinal = 1;

% target set = set of collisions: relative position between quadrotors
x0 = [0; 0; 0; 0];
x0_width = [0.5; 0.5; 0.5; 0.5];
Xend_interval = interval(x0-x0_width,x0+x0_width);
params.R0 = polytope(Xend_interval);

% set of controllable inputs: acceleration in x and y of player 1
params.U = zonotope(interval([-0.5;-0.1],[0.1;0.5]));

% set of uncontrollable disturbances: acceleration in x and y of player 2
params.W = zonotope(interval([0;-0.2],[0.2;0]));

% time step size
options.timeStep = 0.01;

% backward analysis
options.linAlg = 'inner:EA:timeinterval';
R_EA = reachBackward(sys,params,options);
options.linAlg = 'outer:AE:timeinterval';
R_AE = reachBackward(sys,params,options);


% visualization
% projDims = {[1,2],[3,4]};
% for p=1:2
%     subplot(1,2,p); hold on; box on;
%     axis([-2,2,-2,2]);
%     useCORAcolors("CORA:contDynamics",2);
%     plot(R_EA,projDims{p});
%     plot(R_AE,projDims{p});
%     plot(Xend_interval,projDims{p},'k','FaceColor','w');
% end


% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
