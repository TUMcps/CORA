function res = sequentialTrajectories(traj_in, length)
% sequentialTrajectories - Creates sequential trajectories of a given length, 
% starting at every instant of the original trajectory (see Def. 7 in [1])
%
% Syntax:
%    obj = sequentialTrajectories(obj,length)
%
% Inputs:
%    obj - trajectory object
%    length - number of samples for each trajectory
%
% Outputs:
%    res - array of trajectory objects
%
% References:
%    [1] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 2022
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff, Stefan Liu, Laura Luetzow
% Written:       27-August-2025             
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% sanity checks
% check if trajectory object is empty
if isempty(traj_in)
    res = traj_in;
    return
end
% check if length of trajectory is long enough
nrOfTrajectories = traj_in.n_k-length+1;
if nrOfTrajectories < 1
    res = traj_in;
    return
end

% initialize result
res(traj_in.n_s*(traj_in.n_k-length+1),1) = trajectory();
counter = 1;

% loop for each new trajectory
for s = 1:traj_in.n_s
    for i = 0:traj_in.n_k-length
        u = []; x = []; y = []; t = []; a = []; loc = [];

        % update inputs
        if ~isempty(traj_in.u)
            u = traj_in.u(:,i+1:i+length);
        end
        % update states (if state trajectory is provided)
        if ~isempty(traj_in.x) && size(traj_in.x,2) > 1
            x = traj_in.x(:,i+1:i+length,s);
        end
        % update measured outputs
        if ~isempty(traj_in.y)
            y = traj_in.y(:,i+1:i+length,s);
        end
        % update time
        if ~isempty(traj_in.t)
            t = traj_in.t(:,i+1:i+length);
        end
        % update algebraic variables
        if ~isempty(traj_in.a)
            a = traj_in.a(:,i+1:i+length,s);
        end
        % update locations
        if ~isempty(traj_in.loc)
            loc = traj_in.loc(:,i+1:i+length);
        end
        % save result in cell array
        res(counter) = trajectory(u,x,y,t,[],loc,a);
        counter = counter + 1;
    end
end
end

% ------------------------------ END OF CODE ------------------------------
