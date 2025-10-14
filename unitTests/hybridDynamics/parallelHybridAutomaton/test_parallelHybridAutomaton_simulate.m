function res = test_parallelHybridAutomaton_simulate
% test_parallelHybridAutomaton_simulate - test function for simulate
%
% Syntax:
%    res = test_parallelHybridAutomaton_simulate
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

% Authors:       Mark Wetzlinger
% Written:       20-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system dynamics
pHA = roomHeatingParallel();

% simulation
params.startLoc = [1;1];
params.finalLoc = [0;0];
params.tStart = 0;
params.tFinal = 3;
params.x0 = [20.5;20.5];
params.u = 4;

% simulation
[t,x,loc] = simulate(pHA,params);

% must be of same length
assert(size(t,2) == size(x,2) && size(t,2) == size(loc,2));

% check first and final entry in time vector
assert(withinTol(t(1),params.tStart) ...
    && withinTol(t(end),params.tFinal));
% check initial state
assert(compareMatrices(x(:,1),params.x0));
% check start location
assert(compareMatrices(loc(:,1),params.startLoc));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
