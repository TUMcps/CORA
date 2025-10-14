function res = example_linearSysDT_identify
% example_linearSysDT_identify - example for identification of a linear 
%   discrete-time system from trajectory data
%
% Syntax:
%    res = example_linearSysDT_identify
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
% See also: linearSysDT/identify

% Authors:       Niklas Kochdumper
% Written:       17-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Original System ---------------------------------------------------------

A = [0.72 0.36; -0.18 1.08];
dt = 0.2;

sysOrig = linearSysDT(A,dt);


% Trajectory Data ---------------------------------------------------------

params.tFinal = 5;
params.R0 = zonotope(interval([-10;-10],[10;10]));

options.points = 10;

traj = simulateRandom(sysOrig,params,options);
traj = addNoise(traj,0.2,'x','gaussian');


% System Identification ---------------------------------------------------

sys = linearSysDT.identify(traj);


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
