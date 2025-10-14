function res = example_nonlinearSysDT_identify
% example_nonlinearSysDT_identify - example for identification of a 
%   nonlinear discrete-time system from trajectory data
%
% Syntax:
%    res = example_nonlinearSysDT_identify
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] C.V. Rao and et al. "Constrained State Estimation for Nonlinear
%        Discrete-Time Systems: Stability and Moving Horizon 
%        Approximations", IEEE Transactions on Automatic Control 2003
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT/identify

% Authors:       Niklas Kochdumper
% Written:       17-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Original System ---------------------------------------------------------

% system from Equation (9) in [1]
f = @(x,u) [0.99*x(1) + 0.2*x(2); ...
            -0.1*x(1) + 0.5*x(2)/(1+x(2)^2)];
dt = 1;

sysOrig = nonlinearSysDT(f,dt);


% Trajectory Data ---------------------------------------------------------

params.tFinal = 50;
params.R0 = zonotope(interval([-10;-10],[10;10]));

options.points = 10;

traj = simulateRandom(sysOrig,params,options);
traj = addNoise(traj,0.1,'x','gaussian');


% System Identification ---------------------------------------------------

sys = nonlinearSysDT.identify(traj);


% Simulation --------------------------------------------------------------

traj_(length(traj),1) = trajectory();

for i = 1:length(traj)

    simOpts.x0 = traj(i).x(:,1);
    simOpts.tFinal = params.tFinal;

    [t,x] = simulate(sys,simOpts);

    traj_(i) = trajectory([],x,[],t);
end


% Visualization -----------------------------------------------------------

figure; hold on; box on;
plot(traj);
plot(traj_);

res = true;
end

% ------------------------------ END OF CODE ------------------------------
