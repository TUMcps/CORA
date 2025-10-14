function traj = uminus(traj)
% uminus - Overloads the unary '-' operator
%
% Syntax:
%    traj = -traj
%    traj = uminus(traj)
%
% Inputs:
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

traj = -1 * traj;

end

% ------------------------------ END OF CODE ------------------------------
