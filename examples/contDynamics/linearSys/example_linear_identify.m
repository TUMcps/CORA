function res = example_linear_identify
% example_linear_identify - example for identification of a linear system
%   from trajectory data
%
% Syntax:
%    res = example_linear_identify
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

% Authors:       Niklas Kochdumper
% Written:       06-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Original System ---------------------------------------------------------

A = [-0.7 -2; 2 -0.7];

sysOrig = linearSys(A);


% Trajectory Data ---------------------------------------------------------

params.tFinal = 5;
params.R0 = zonotope(interval([-10;-10],[10;10]));

options.points = 10;

traj = simulateRandom(sysOrig,params,options);
traj = addNoise(traj,0.1,'x','gaussian');


% System Identification ---------------------------------------------------

sys = linearSys.identify(traj);


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
