function res = example_linear_reachBackward_quadrotor6D
% example_linear_reachBackward_quadrotor6D - example for backward
%    reachability analysis
%
% Syntax:
%    res = example_linear_reachBackward_quadrotor6D
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] I.M. Mitchell, J. Budzis, A. Bolyachevets. "Invariant, Viability
%        and Discriminating Kernel Under-Approximation via Zonotope
%        Scaling". arXiv, 2019.
%    [2] A. Sasfi, M. Zeilinger, J. Koehler. "Robust adaptive MPC using
%        control contraction metrics", Automatica, 2023.
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

% linearization point: only roll and total thrust required
xlin = [0; 0; 0; 0; 0; 0];
ulin = [9.81/(0.89/1.4); 0];

% init model
[A,B,E] = quadrotor6D(xlin,ulin);
sys = linearSys(A,B,[],[],[],[],E);

% target set for ground collision avoidance:
% - horizontal position/velocity does not matter
% - vertical position not below zero (+ safe distance)
% - horizontal velocity not below zero
% - roll/roll velocity from state constraints [1, (44)]
H = [1 -1 0  0 1  -1 0  0 0  0  0  0 0  0  0;
     0  0 1 -1 -2 -2 0  0 0  0  10 0 0  0  0;
     0  0 0  0 0   0 1 -1 0  0  0  0 0  0  0;
     0  0 0  0 0   0 0  0 1 -1 -1  0 0  0  0;
     0  0 0  0 0   0 0  0 0  0  0  1 -1 0  0;
     0  0 0  0 0   0 0  0 0  0  0  0 0  1 -1]';
d = [0.5 0.5 0.1 0 0.3 0.3 1 1 0 1 1 pi/15 pi/15 pi/2 pi/2]';
params.R0 = polytope(H,d);

params.tStart = 0;
params.tFinal = 0.5;

% algorithm and time step size
options.linAlg = 'outer:AE:timeinterval';
options.timeStep = 0.0025;

% scaling factors for different input/disturbance capacities
scaleU = [1; 1; 2];
lbw = {[-10*0.2760;-0.3668],[-0.2760;-0.3668],[-0.2760;-0.3668]};
ubw = {[10*0.2760;0.3668],[0.2760;0.3668],[0.2760;0.3668]};

R = cell(3,1);
for s=1:3
    % set of controllable inputs [1, (45)]: total thrust and desired roll angle
    params.U = zonotope([ulin(1);xlin(5)],diag([scaleU(s)*1.5; pi/6]));
    
    % set of uncontrollable disturbance [2, Sec. 4]: wind
    params.W = zonotope(interval(lbw{s},ubw{s}));

    % backward reachability analysis
    R{s} = reachBackward(sys,params,options);
end


% visualization
% figure;
% projDims = {[1,3],[2,4],[1,2]};
% for p=2:length(projDims)
%     subplot(1,3,p); hold on; box on;
%     useCORAcolors("CORA:contDynamics",3);
%     for s=1:3
%         plot(R,projDims{p});
%     end
% end


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
