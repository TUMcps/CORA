function traj = simResult(varargin)
% simResult - (deprecated) object constructor for simResult
%
% Syntax:
%    traj = simResult(x,t)
%    traj = simResult(x,t,loc)
%    traj = simResult(x,t,{},y)
%    traj = simResult(x,t,{},y,a)
%
% Inputs:
%    x - cell-array storing the states of the simulated trajectories, where
%        each trajectory is a matrix of dimension [N,n]
%    t - cell-array storing the time points for the simulated trajectories
%        as vectors of dimension [N,1]
%    loc - double-array storing the locations for the simulated trajectories
%    y - cell-array storing the output of the simulated trajectories, where
%        each trajectory is a matrix of dimension [N,o]
%    a - cell-array storing the algebraic variables of the simulated
%        trajectories (only for nonlinDASys), where each trajectory is a
%        matrix of dimension [N,p]
%    (where N ... number of simulated trajectories, n ... state dimension,
%     o ... output dimension, p ... number of algebraic variables)
%
% Outputs:
%    traj - trajectory object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: trajectory

% Authors:       Laura Luetzow
% Written:       25-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

CORAwarning("CORA:deprecated", "class", "simResult", ...
    "CORA v2026", ...
    "Please use trajectory instead. The object is returned as a trajectory.",...
    "The new trajectory class unifies the classes simResult and testCase.")

% avoid empty instantiation
if nargin == 0
    % empty object
    return
end
assertNarginConstructor(1:5,nargin);

if nargin == 1 && isa(varargin{1},'simResult')
    % convert to trajectory
    traj = simResult2traj(varargin{1});
    return
end
x = {}; t = {}; loc = {}; y = {}; a = {};
[x,t,loc,y,a] = setDefaultValues({x,t,loc,y,a},varargin);

% create dummy simResult object and convert to trajectory
simRes = struct('x', {x}, 't', {t}, 'loc', {loc}, 'y', {y}, 'a', {a});

traj = simResult2traj(simRes);

end

% ------------------------------ END OF CODE ------------------------------
