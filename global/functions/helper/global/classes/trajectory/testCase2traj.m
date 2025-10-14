function traj = testCase2traj(testC)
% testCase2traj - function to convert testCase object to trajectory object
%
% Syntax:
%    traj = testCase2traj(testC)
%
% Inputs:
%    testC - testCase object (deprecated) or dummy struct
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
% Written:       26-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


dt = testC.sampleTime;
if (size(testC.u,3) == 1 || all(diff(testC.u, [], 3) == 0, 'all')) && ...
        (size(testC.initialState,3) == 1 || all(diff(testC.initialState, [], 3) == 0, 'all'))
    % return one trajectory with multiple realizations starting from the
    % same initial state and subject to the same inputs
    u = testC.u(:,:,1)';
    if isempty(testC.x)
        % if no sate trajectory measured, use initial state
        x = testC.initialState(:,:,1);
    else
        x = permute(testC.x, [2 1 3]);
    end
    y = permute(testC.y, [2 1 3]);
    if isprop(testC, 'model')
        traj = trajectory(u,x,y,[],dt,[],[],testC.name, testC.model);
    else
        traj = trajectory(u,x,y,[],dt,[],[],testC.name);
    end

else
    % return multiple trajectories

    n_m = max([size(testC.x,3), size(testC.y,3)]);
    traj(n_m,1) = trajectory();
    for m = 1:n_m
        % convert inputs
        if size(testC.u,3) > 1
            u = testC.u(:,:,m)';
        else
            u = testC.u';
        end

        % convert states
        if size(testC.x,3) > 1
            x = testC.x(:,:,m)';
        else
            x = testC.x';
        end
        if isempty(x)
            % if no sate trajectory measured, use initial state
            if size(testC.initialState,3) > 1
                x = testC.initialState(:,:,m);
            else
                x = testC.initialState;
            end
        end

        % convert outputs
        if size(testC.y,3) > 1
            y = testC.y(:,:,m)';
        else
            y = testC.y';
        end

        % create trajectory object
        if isprop(testC, 'model')
            traj(m) = trajectory(u,x,y,[],dt,[],[], testC.name, testC.model);
        else
            traj(m) = trajectory(u,x,y,[],dt,[],[], testC.name);
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
