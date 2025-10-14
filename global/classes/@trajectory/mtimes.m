function traj = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the states of the simulated trajectories
%
% Syntax:
%    traj = mtimes(A,traj)
%
% Inputs:
%    A - numeric matrix
%    traj - trajectory object
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

[traj,A] = findClassArg(factor1,factor2,'trajectory');

if ~isnumeric(A)
    throw(CORAerror("CORA:noops",A,traj))
end

% compute linear map for each trajectory
for r=1:length(traj)
    traj(r).x = pagemtimes(A, traj(r).x);
end

end

% ------------------------------ END OF CODE ------------------------------
