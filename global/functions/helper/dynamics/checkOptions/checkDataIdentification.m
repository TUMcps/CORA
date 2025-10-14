function [traj,options] = checkDataIdentification(varargin)
% checkDataIdentification - Checks if the trajectory data provided for
%   system identification has the correct format
%
% Syntax:
%    [traj,options] = checkDataIdentification(traj)
%    [traj,options] = checkDataIdentification(traj,options)
%    [traj,options] = checkDataIdentification(x,t)
%    [traj,options] = checkDataIdentification(x,t,options)
%    [traj,options] = checkDataIdentification(x,t,u)
%    [traj,options] = checkDataIdentification(x,t,u,options)
%
% Inputs:
%    traj - object of class "trajectory" storing the trajectory data
%    x - cell-array storing the states of the simulated trajectories, where
%        each trajectory is a matrix of dimension [n,N]
%    t - cell-array storing the time points for the simulated trajectories
%        as vectors of dimension [1,N]
%    u - cell-array storing the inputs for the simulated trajectories
%        as vectors of dimension [m,N-1] or [m,N]
%    options - algorithm options for system identification
%       .alg: algorithm that is used ('dmd': Dynamic Mode Decomposition, 
%             'opt': optimization)
%
% Outputs:
%    traj - object of class "trajectory" storing the data
%    options - algorithm options for system identification
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       05-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    options = [];

    % check if number of input arguments is correct
    try
        narginchk(1,4);
    catch ex
        throwAsCaller(ex)
    end

    % check which format the input is in
    if isa(varargin{1},'trajectory')

        % check if the number of input arguments is correct
        try
            narginchk(1,2);
        catch ex
            throwAsCaller(ex)
        end

        if nargin > 1
            options = varargin{2}; posOptions = 'second';
        end

    else
        
        % check if the number of input arguments is correct
        try
            narginchk(2,4);
        catch ex
            throwAsCaller(ex)
        end

        % go over the different numbers of input arguments
        x = varargin{1}; t = varargin{2}; u = [];

        if nargin == 4
            u = varargin{3}; options = varargin{4}; posOptions = 'forth'; 
        elseif nargin == 3

            if isstruct(varargin{3})
                options = varargin{3}; posOptions = 'third';
            else
                u = varargin{3};
            end
        end

        % convert to cell-arrays
        if ~iscell(x)
            x = {x};
        end

        if ~iscell(t)
            t = {t};
        end

        if ~isempty(u) && ~iscell(u)
            u = {u};
        end
    end

    % check if the inputs are correct
    if ~isa(varargin{1},'trajectory')
        aux_checkTrajectoryFormat(x,t,u);
    end

    % check if options are correct
    if ~isempty(options) && ~isstruct(options)
        throwAsCaller(CORAerror('CORA:wrongValue',posOptions,'struct'));
    end

    % bring data to a common format
    if isa(varargin{1},'trajectory')

        traj_old = varargin{1};
        traj(sum([traj_old.n_s]),1) = trajectory();
        counter = 1;
        for i = 1:length(traj_old)
            for s=1:traj_old(i).n_s
                % input
                u = [];
                if ~isempty(traj_old(i).u)
                    u = traj_old(i).u;
                end
                % output
                y = [];
                if ~isempty(traj_old(i).y)
                    y = traj_old(i).y(:,:,s);
                end
                % state
                x = traj_old(i).x(:,:,s);
                % time
                t = traj_old(i).t;
                %create trajectory
                traj(counter) = trajectory(u,x,y,t);
                counter = counter + 1;
            end
        end

    else
        % signals are provided separately
        traj(length(x),1) = trajectory();
        for i = 1:length(x)
            x_i = x{i};
            t_i = t{i};
            u_i = [];
            if ~isempty(u)
                u_i = u{i};
            end
            traj(i) = trajectory(u_i,x_i,[],t_i);
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function aux_checkTrajectoryFormat(x,t,u)

    % check correct format of simulated trajectories
    if CHECKS_ENABLED

        % - time has to be non-empty, the entries have to be nonnegative
        %   and monotonically increasing
        if (isempty(t) && ~isempty(x)) || ~iscell(t) || ...
                any(cellfun(@(x_) size(x_,1),t,'UniformOutput',true) ~= 1) || ...
                any(cellfun(@(x_) any(x_ < 0),t,'UniformOutput',true)) || ...
                any(cellfun(@(x_) any(diff(x_ < 0)),t,'UniformOutput',true))
            throwAsCaller(CORAerror('CORA:wrongValue','second', ...
                'non-negative and monotonically increasing vector or cell-array of vectors'));
        end

        tVeclength = cellfun(@(x_) length(x_),t,'UniformOutput',true);
    
        % number of steps of state vector has to be equal length of time vector
        proplength = cellfun(@(x_) size(x_,2),x,'UniformOutput',true);
        if ~isempty(proplength) && any(proplength ~= tVeclength)
                throwAsCaller(CORAerror('CORA:wrongValue','first', ...
                    'The length of each run has to match the length of the time vector.'));
        end

        % number of steps of input vector has to be equal length of time vector
        if ~isempty(u)
            proplength = cellfun(@(x_) size(x_,2),u,'UniformOutput',true);
            if ~isempty(proplength) && any(proplength ~= tVeclength) && any(proplength+1 ~= tVeclength)
                    throwAsCaller(CORAerror('CORA:wrongValue','third', ...
                        'The length of each run has to match the length of the time vector.'));
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
