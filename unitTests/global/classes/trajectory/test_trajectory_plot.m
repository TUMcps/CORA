function res = test_trajectory_plot
% test_trajectory_plot - unit test function for plot
%
% Syntax:
%    res = test_trajectory_plot()
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

% Authors:       Mark Wetzlinger, Laura Luetzow
% Written:       07-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init result
res = true;
figure;


% simulations from linear discrete-time system
A = [0.9810    0.0143    0.0262   -0.0140;
    0.0079    1.0133    0.0267    0.0416;
    0.0029   -0.0054    0.9680    0.0188;
    0.0503   -0.0242    0.0411    0.9877];
B = [-0.0246    0.0003    0.0372;
    0.0251    0.0121    0.0006;
   -0.0001   -0.0003    0.0125;
   -0.0009   -0.0250    0.0136];
dt = 0.05;
sys = linearSysDT(A,B,dt);
params.R0 = zonotope([10;-5;8;-12],0.5*eye(4));
params.tFinal = 5;


% single point
simOpt.points = 1;
traj = simulateRandom(sys,params,simOpt);

% plot time step
plot(traj,[2,4]);
plot(traj,[1,3],'r');

% multiple points
simOpt.points = 10;
traj = simulateRandom(sys,params,simOpt);

% plot time step
plot(traj);
plot(traj,[1,3],'r');

% simulations from nonlinear continuous-time system
f = @(x,u) [-x(1)^2*x(2) - u(1); -exp(x(2))];
sys = nonlinearSys(f,2,1);
params.R0 = zonotope([10;5],0.05*eye(2));
params.U = zonotope(0,0.02);
params.tFinal = 2;


% single point
simOpt.points = 1;
traj = simulateRandom(sys,params,simOpt);

% plot time step
plot(traj);
plot(traj,[1,2],'r');

% multiple points
simOpt.points = 10;
traj = simulateRandom(sys,params,simOpt);

% plot time step
plot(traj,[1,2]);
plot(traj,[1,2],'r');

% close figure
close

% ------------------------------ END OF CODE ------------------------------
