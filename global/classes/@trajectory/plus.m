function traj = plus(traj,v)
% plus - Overloaded '+' operator for the states of the simulated trajectories
%
% Syntax:
%    traj = plus(traj,v)
%
% Inputs:
%    traj - trajectory object
%    S - contSet object or numerical vector
%
% Outputs:
%    traj - transformed trajectory object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Laura Luetzow
% Written:       07-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get reachSet object
[traj,v] = findClassArg(traj,v,'trajectory');

if ~isnumeric(v)
    throw(CORAerror("CORA:noops",traj,v))
end

% add to all state trajectories
for r=1:length(traj)
    traj(r).x = traj(r).x + v;
end

end

% ------------------------------ END OF CODE ------------------------------
