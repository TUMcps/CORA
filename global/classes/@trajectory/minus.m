function traj = minus(traj1,traj2)
% minus - Overloaded '-' operator for the states of the simulated trajectories
%
% Syntax:
%    traj = minus(traj,factor2)
%
% Inputs:
%    traj1 - numeric or trajectory object
%    factor2 - numeric or trajectory object 
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

traj = traj1 + (-traj2);

end

% ------------------------------ END OF CODE ------------------------------
