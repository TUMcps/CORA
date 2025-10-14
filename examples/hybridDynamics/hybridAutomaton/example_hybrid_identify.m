function res = example_hybrid_identify
% example_hybrid_identify - example for identification of a hybrid
%   automaton from trajectory data using the method in [1]
%
% Syntax:
%    res = example_hybrid_identify
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%   [1] N. Kochdumper and et al. "Robust Identification of Hybrid Automata 
%       from Noisy Data", HSCC 2025
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

rng(1)
HAorig = bouncing_ball(-0.75);


% Trajectory Data ---------------------------------------------------------

% simulation parameter
params.R0 = interval([1;-1],[2;1]);
params.startLoc = 1;
params.tFinal = 2;

% simulate benchmark
traj = simulateRandom(HAorig,params,struct('points',10)); 

% add noise to the data
traj = addNoise(traj,0.1,'x','gaussian');


% System Identification ---------------------------------------------------

HA = hybridAutomaton.identify(traj);


% Simulation --------------------------------------------------------------

traj_(length(traj),1) = trajectory();

for i = 1:length(traj)

    simOpts.x0 = traj(i).x(:,1);
    simOpts.tFinal = params.tFinal;
    simOpts.startLoc = bestInitialMode(HA,traj(i).x,traj(i).t);

    [t_,x_] = simulate(HA,simOpts);

    traj_(i) = trajectory([],x_,[],t_);
end


% Visualization -----------------------------------------------------------

figure; hold on; box on;
plotOverTime(traj(end),1);
plotOverTime(traj_(end),1);

res = true;

end

% ------------------------------ END OF CODE ------------------------------
