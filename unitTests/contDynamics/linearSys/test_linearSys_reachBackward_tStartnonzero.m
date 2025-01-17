function res = test_linearSys_reachBackward_tStartnonzero
% test_linearSys_reachBackward_tStartnonzero - unit test for backward
%    reachability analysis (time-interval solution), where the start time
%    is non-zero
%
% Syntax:
%    res = test_linearSys_reachBackward_tStartnonzero
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
% Written:       29-November-2024
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

% target set = set of collisions: relative position between quadrotors
params.R0 = polytope(interval(-ones(4,1),ones(4,1)));
% set of controllable inputs: acceleration in x and y of player 1
params.U = zonotope(zeros(2,1),eye(2));
% set of uncontrollable disturbances: acceleration in x and y of player 2
params.W = zonotope(zeros(2,1),eye(2));

% time step size divides [0,tStart] and [tStart,tFinal] in integer steps
params.tStart = 0.20;
params.tFinal = 0.30;
options.timeStep = 0.05;

% EA backward reachability
options.linAlg = 'inner:EA:timeinterval';
R = reachBackward(sys,params,options);

% AE backward reachability
options.linAlg = 'outer:AE:timeinterval';
R = reachBackward(sys,params,options);


% time step size does not divide [0,tStart] in integer steps
params.tStart = 0.10;
params.tFinal = 0.18;
options.timeStep = 0.04;

% EA backward reachability
options.linAlg = 'inner:EA:timeinterval';
R = reachBackward(sys,params,options);

% AE backward reachability
options.linAlg = 'outer:AE:timeinterval';
R = reachBackward(sys,params,options);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
