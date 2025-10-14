function traj = simResult2traj(simRes)
% simResult2traj - function to convert simResult object to trajectory object
%
% Syntax:
%    traj = simResult2traj(simRes)
%
% Inputs:
%    simRes - simResult object (deprecated) or dummy struct
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
% Written:       25-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n_m = max([length(simRes.x), length(simRes.y), length(simRes.t), length(simRes.a)]);
traj(n_m,1) = trajectory();

for m = 1:n_m
    u = [];
    % convert states
    if ~isempty(simRes.x)
        x = simRes.x{m}';
    else
        x = [];
    end

    % convert outputs
    if ~isempty(simRes.y)
        y = simRes.y{m}';
    else
        y = [];
    end

    % convert time vector
    if ~isempty(simRes.t)
        t = simRes.t{m}';
    else
        t = [];
    end

    % convert algebraic variables
    if ~isempty(simRes.a)
        a = simRes.a{m}';
    else
        a = [];
    end

    % convert locations
    if ~isempty(simRes.loc)
        n_k = max([size(x,2), size(y,2), size(t,2), size(a,2)]);
        loc = repmat(simRes.loc(m), 1, n_k);
    else
        loc = [];
    end

    % create trajectory object
    traj(m) = trajectory(u,x,y,t,[],loc,a);
end
end

% ------------------------------ END OF CODE ------------------------------
