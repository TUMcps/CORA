function res = isemptyobject(traj)
% isemptyobject - checks if a trajectory object is empty
%
% Syntax:
%    res = isemptyobject(traj)
%
% Inputs:
%    traj - trajectory object
%
% Outputs:
%    res - true/false
%
% Example:
%    traj = trajectory();
%    isemptyobject(traj)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger, Laura Luetzow
% Written:       07-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[r,c] = size(traj);
res = true;

% loop over all objects
for i=1:r
    for j=1:c
        % check time
        res = res && isempty(traj(i,j).t) && isempty(traj(i,j).u) ...
            && isempty(traj(i,j).x) && isempty(traj(i,j).y) ...
            && isempty(traj(i,j).loc) && isempty(traj(i,j).a);
    end
end

% ------------------------------ END OF CODE ------------------------------
