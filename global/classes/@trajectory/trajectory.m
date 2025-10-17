classdef trajectory
% trajectory - class that stores trajectory data
%
% Syntax:
%    traj = trajectory(u,x,y,t)
%    traj = trajectory(u,x,y,t,dt)
%    traj = trajectory(u,x,y,t,dt,loc)
%    traj = trajectory(u,x,y,t,dt,[])
%    traj = trajectory(u,x,y,t,dt,[],a) 
%    traj = trajectory(u,x,y,t,dt,[],a,name)
%    traj = trajectory(u,x,y,t,dt,[],a,name,model)
%
% Inputs:
%    u - (dim_u x n_k) array of inputs
%    x - (dim_x x n_k x n_s) array of states
%    y - (dim_y x n_k x n_s) array of the measured outputs
%    t - (1     x n_k) array of time points 
%    dt - sampling time
%    a - (dim_a x n_k x n_s) array of algebraic variables
%        (only for nonlinDASys)
%    loc - (dim_l x n_k) array storing the locations for the trajectories
%    name - name of the test case
%    (where n_k ... number of time steps, n_s ... number of realizations 
%     starting from the same initial state and subject to the same inputs), 
%     dim_u ... input dimension, dim_x ... state dimension,
%     dim_y ... output dimension, dim_a ... algebraic variable dimension,
%     dim_l ... location dimension)
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet, simulateRandom

% Authors:       Laura Luetzow
% Written:       25-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    
    u = [];             % inputs of the trajectories
    x = [];             % states of the trajectories
    y = [];             % outputs of the trajectories
    t = [];             % time points of the trajectories
    dt = 0;             % average time step inbetween time points
    a = [];             % algebraic variables of the trajectories
    loc = [];            % index of the locations (0 for contDynamics)
    name = "";          % name of the trajectories
    constant_dt = false;% constant sampling time
    n_k = 0;            % number of time steps of trajectory
    n_s = 0;            % number of trajectories saved
    model               % model that created the trajectories
    
end

methods
    
    % class constructor
    function traj = trajectory(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            % empty object
            return
        end
        assertNarginConstructor(1:9,nargin);

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'trajectory')
            % copy constructor
            traj = varargin{1};
        end
        if nargin == 9
            % system model is provided
            traj.model = varargin{9};
            varargin = varargin(1:8);
        end
        
        % 2. parse input arguments: varargin -> vars
        [u,x,y,t,dt,loc,a,name] = aux_parseInputArgs(varargin{:});

        n_k = max([size(u,2), size(x,2), size(y,2), size(t,2), size(a,2)]);
        n_s = max([size(x,3), size(y,3), size(a,3)]);
        if isempty(t) && ~isempty(dt)
            % construct time points from sampling time
            t = 0:dt:(n_k-1)*dt;
        elseif ~isempty(t) && isempty(dt) && length(t) > 1
            % compute average sampling time from time points
            dt = mean(diff(t,[],2),'all');
        end

        % 3. check correctness of input arguments
        aux_checkInputArgs(u,x,y,t,dt,loc,a,name,n_k,n_s,nargin);

        % 4. assign properties
        traj.x = x;
        traj.y = y;
        traj.u = u;
        traj.t = t;
        traj.dt = dt;
        traj.loc = loc;
        traj.a = a;
        traj.name = name;
        traj.constant_dt = all(withinTol(diff(t,[],2), dt, 1e-6), 'all'); % TO-DO: hybrides System?
        traj.n_k = n_k;
        traj.n_s = n_s;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [u,x,y,t,dt,loc,a,name] = aux_parseInputArgs(varargin)

    u = []; x = []; y = []; t = []; dt = []; loc = []; a = []; name = "";
    [u,x,y,t,dt,loc,a,name] = setDefaultValues({u,x,y,t,dt,loc,a,name},varargin);
end

function aux_checkInputArgs(u,x,y,t,dt,loc,a,name,n_k,n_s,n_in)

% check correct format of trajectories
if CHECKS_ENABLED && n_in > 0

    % input must be empty or have dimensions
    % (? x n_k x 1)
    if ~isempty(u) && ((size(u,2) ~= n_k && size(u,2) ~= n_k-1) || size(u,3) ~= 1)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Property .u has the wrong format.'));
    end

    % state must be empty or have dimensions
    % (? x n_k x n_s) with the same initial states or 
    % (? x 1   x 1) if only initial state is known 
    if ~isempty(x) && ~((size(x,2) == n_k && size(x,3) == n_s && ...
            all(abs(diff(x(:, 1, :), [], 3)) < 1e-6 , 'all')) || ...
            (size(x,2) == 1 && size(x,3) == 1))
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Property .x has the wrong format.'));
    end

    % output dimensions must be empty or 
    % (? x n_k x n_s)
    if ~isempty(y) && (size(y,2) ~= n_k || size(y,3) ~= n_s)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Property .y has the wrong format.'));
    end
    
    % time points have to be nonnegative and monotonically increasing, and
    % must be empty or have dimensions
    % (1 x n_k x 1) 
    if ~isempty(t) && (size(t,1) ~= 1 || size(t,2) ~= n_k || ...
            size(t,3) ~= 1 || any(t < 0,'all') || any(diff(t,[],2) < 0,'all'))
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Property .t has the wrong format.'));
    end

    % time points and must match the average sampling time
    if ~isempty(t) && ~isempty(dt) && n_k > 1 && ~withinTol(dt, mean(diff(t,[],2),'all'), 1e-6)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Time points and average sampling time does not match.'));
    end

    % locations must be empty or have dimensions
    % (? x n_k x 1) 
    if ~isempty(loc) && (size(loc,2) ~= n_k || size(loc,3) ~= 1)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Property .loc has the wrong format.'));
    end

    % algebraic variables must be [] or have dimensions
    % (? x n_k x n_s)
    if ~isempty(a) && (size(a,2) ~= n_k || size(a,3) ~= n_s)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Property .a has the wrong format.'));
    end
end

end

% ------------------------------ END OF CODE ------------------------------
