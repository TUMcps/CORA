function res = example_nonlinear_identify
% example_nonlinear_identify - example for identification of a nonlinear 
%   system from trajectory data
%
% Syntax:
%    res = example_nonlinear_identify
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
% Written:       13-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Original System ---------------------------------------------------------

f = @(x,u) [x(2); ...
            (1-x(1)^2)*x(2)-x(1)];

sysOrig = nonlinearSys(f);


% Trajectory Data ---------------------------------------------------------

params.tFinal = 10;
params.R0 = zonotope(interval([-3;-3],[3;3]));

options.points = 10;

traj = simulateRandom(sysOrig,params,options);


% System Identification ---------------------------------------------------

sys = nonlinearSys.identify(traj);


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
