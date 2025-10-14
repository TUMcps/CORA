function res = test_simResult
% test_simResult - unit test function for (deprecated) simResult 
%   constructor returning a trajectory object
%
% Syntax:
%    res = test_simResult()
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

% Authors:       Laura Luetzow
% Written:       26-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize some sets for instantiations
steps = 10;
runs = 5;

nrStates = 4; nrOutputs = 2; nrAlgebraic = 3;
x = cell(runs,1); 
y = cell(runs,1); 
a = cell(runs,1); 
t = cell(runs,1);
loc = ones(runs,1);

% instantiate trajectories
for i=1:runs
    x{i,1} = randn(steps+1,nrStates);
    y{i,1} = randn(steps+1,nrOutputs);
    a{i,1} = randn(steps+1,nrAlgebraic);
    t{i,1} = cumsum([0;rand(steps,1)]);
end

% correct instantiations
traj = simResult(x,t);
assert(isa(traj, 'trajectory'))
traj = simResult(x,t,loc);
assert(isa(traj, 'trajectory'))
traj = simResult(x,t,{},y);
assert(isa(traj, 'trajectory'))
traj = simResult(x,t,{},y,a);
assert(isa(traj, 'trajectory'))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
