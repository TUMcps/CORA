function traj = uniformTimeStepSize(traj,dt)
% uniformTimeStepSize - convert the trajectory data to a uniform time step 
%   size
%
% Syntax:
%    traj = uniformTimeStepSize(traj, dt)
%
% Inputs:
%    traj - trajectory object
%    dt - time step size
%
% Outputs:
%    traj - new trajectory object with time step size dt
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: trajectory

% Authors:       Niklas Kochdumper, Laura Luetzow
% Written:       28-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% select time step size
if nargin < 2
    t = [];

    for i = 1:length(traj)
        t = [t diff(traj(i).t)];
    end

    dt = mean(t);
end

% convert all measured traces to discrete time
for i = 1:length(traj)

    d = diff(traj(i).t);
    if all(abs(d - dt) < eps)
        continue;
    end

    t = (0:floor(traj(i).t(end)/dt))*dt;
    [~,ind] = unique(traj(i).t);

    % interpolate states
    x = [];
    if ~isempty(traj(i).x)
        x = interp1(traj(i).t(ind),traj(i).x(:,ind)',t,'linear','extrap')';
        if size(x,1) ~= size(traj(i).x,1)
            % interp1 returns (number of steps x 1) vector if number of states = 1,
            % otherwise (number of states x number of steps) array 
            x = x';
        end
    end

    % interpolate inputs
    u = [];
    if ~isempty(traj(i).u)
        if size(traj(i).u,2) < size(traj(i).x,2)
            traj(i).u = [traj(i).u traj(i).u(:,end)];
        end
        u = interp1(traj(i).t(ind),traj(i).u(:,ind)',t,'linear','extrap')';
        if size(u,1) ~= size(traj(i).u,1)
            u = u';
        end
    end

    % interpolate outputs
    y = [];
    if ~isempty(traj(i).y)
        y = interp1(traj(i).t(ind),traj(i).y(:,ind)',t,'linear','extrap')';
        if size(y,1) ~= size(traj(i).y,1)
            y = y';
        end
    end
    traj(i) = trajectory(u,x,y,t,dt);
end

% ------------------------------ END OF CODE ------------------------------
