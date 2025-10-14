function y_a = computeOutputDev(traj,sys)
% computeOutputDev - compute measurement deviation y_a by subtracting the 
%   nominal solution from the measurement trajectory
%
% Syntax:
%    y_a = computeOutputDev(traj, sys)
%
% Inputs:
%    traj - trajectory object
%    sys - contDynamics object for computing nominal output trajectory
%
% Outputs:
%    y_a - deviaton between measured outputs and nominal outputs
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: trajectory

% Authors:       Laura Luetzow
% Written:       28-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

y_a = [];
for i = 1:length(traj)
    params.tFinal = sys.dt * size(traj(i).y,2) - sys.dt;
    params.x0 = traj(i).x(:,1,1);
    params.u = traj(i).u;

    % compute nominal outputs via simulation of the dynamics
    [~,~,~,y_nom] = simulate(sys, params);
    y_ai = traj(i).y - y_nom;

    % fill y_a with NaN if concatenated arrays have not the same size
    if size(y_a,2) > size(y_ai,2)
        y_ai = [y_ai NaN(size(y_ai,1), size(y_a,2) - size(y_ai,2), size(y_ai,3))];
    elseif size(y_a,2) < size(y_ai,2) && i > 1
        y_a = [y_a NaN(size(y_a,1), size(y_ai,2) - size(y_a,2), size(y_a,3))];
    end
    y_a = cat(3,y_a,y_ai);
end
end

% ------------------------------ END OF CODE ------------------------------
