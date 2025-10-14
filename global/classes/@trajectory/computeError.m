function err = computeError(traj,sys)
% computeError - compute the error between the system approximation 
%       and the real data
%
% Syntax:
%    y_a = computeError(traj, sys)
%
% Inputs:
%    traj - trajectory object containing the real data
%    sys - system object
%
% Outputs:
%    err - mean squared error
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

err = 0;

% loop over all trajectories
for i = 1:length(traj)

    % simulation options
    simOpts = [];
    simOpts.tStart = traj(i).t(1);
    simOpts.tFinal = traj(i).t(end);

    if isa(sys,'linearARX')
        simOpts.y0 = traj(i).x(:,1:sys.n_p);
    else
        simOpts.x0 = traj(i).x(:,1);
    end
    
    if ~isempty(traj(i).u)
        simOpts.u = traj(i).u;
        if isa(sys,'nonlinearSysDT')
            simOpts.u = simOpts.u(:,1:length(traj(i).t)-1);
        end
    end

    % simulate the system
    if isa(sys,'linearARX')
        [t,~,~,x] = simulate(sys,simOpts);
    else
        [t,x] = simulate(sys,simOpts);
    end

    % compute the error for the current trajectory
    [~,ind] = unique(t);
    x = interp1(t(ind),x(:,ind)',traj(i).t,'linear','extrap')';
    if size(x,1) ~= size(traj(i).x,1)
        % interp1 returns (number of steps x 1) vector if number of states = 1,
        % otherwise (number of states x number of steps) array
        x = x';
    end

    err = err + mean(sum((x-traj(i).x).^2,1));
end
end

% ------------------------------ END OF CODE ------------------------------
