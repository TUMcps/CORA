function sys = identify(varargin)
% identify - Identifies a linear discrete-time system from trajectory data
%
% Syntax:
%    sys = linearSysDT.identify(traj)
%    sys = linearSysDT.identify(traj,options)
%    sys = linearSysDT.identify(x,t)
%    sys = linearSysDT.identify(x,t,options)
%    sys = linearSysDT.identify(x,t,u)
%    sys = linearSysDT.identify(x,t,u,options)
%
% Inputs:
%    traj - object of class "trajectory" storing the trajectory data
%    x - cell-array storing the states of the simulated trajectories, where
%        each trajectory is a matrix of dimension [n,N]
%    t - cell-array storing the time points for the simulated trajectories
%        as vectors of dimension [1,N]
%    u - cell-array storing the inputs for the simulated trajectories
%        as vectors of dimension [m,N-1]
%    options - algorithm options for system identification
%
%       -.dt:   sampling time for the discrete-time system. The default
%               value is the average time step size from the data. 
%       -.alg:  algorithm that is used ('dmd': Dynamic Mode 
%               Decomposition [1], 'opt': optimization). The default value
%               is 'dmd'.
%
% Outputs:
%    sys - identified linearSysDT object
%
% Example: 
%    A = [0 1;-0.5 1];
%    dt = 1;
%    sysOrig = linearSysDT(A,dt);
%
%    simOpts.x0 = [1;-5];
%    simOpts.tFinal = 15;
%    [t,x] = simulate(sysOrig,simOpts);
%
%    sys = linearSysDT.identify(x,t);
%
%    [t_,x_] = simulate(sys,simOpts);
%
%    figure; hold on; box on;
%    plot(x(1,:),x(2,:));
%    plot(x_(1,:),x_(2,:),'LineStyle','--');
%
% References:
%    [1] J.L. Proctor and et al. "Dynamic mode decomposition with control" 
%         SIAM Journal on Applied Dynamical Systems 15, 1 (2016), 142–161
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/identify

% Authors:       Niklas Kochdumper
% Written:       17-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse the input arguments and bring trajectory data to suitable form
    [traj,options] = checkDataIdentification(varargin{:});

    % check the algorithm options
    options = aux_checkOptions(options,traj);

    % convert data to uniform time step size
    traj = uniformTimeStepSize(traj,options.dt);

    % call the selected algorithm
    switch options.alg

        case 'dmd'
            sys = priv_identifyDMD(traj);

        case 'opt'
            sys = priv_identifyOpt(traj);
    end
end


% Auxiliary functions -----------------------------------------------------

function options = aux_checkOptions(options,data)
% check the algorithm settings provided by the user

    % default values
    if ~isfield(options,'alg')
        options.alg = 'dmd';
    end

    if ~isfield(options,'dt')
        options.dt = aux_averageTimeStepSize(data);
    end

    % check user input
    if options.dt <= 0 || ~isscalar(options.dt)
        throw(CORAerror('CORA:wrongFieldValue','options.dt','double > 0'));
    end

    if ~ismember(options.alg,{'dmd','opt'})
        throw(CORAerror('CORA:wrongFieldValue','options.alg', ...
                                                        {'dmd','opt'}));
    end

    % check if there are any redundant options specified
    redundantOptions(options,{'alg','dt'});
end

function dt = aux_averageTimeStepSize(traj)
% compute the average time step size from the trajectory data
    
    % loop over all trajectories
    dt = 0; N = 0;

    for i = 1:length(traj)

        % remove duplicate times
        [~,ind] = unique(traj(i).t);

        % add time steps for current trajectory
        dt = dt + sum(diff(traj(i).t(ind)));
        N = N + length(ind) - 1;
    end

    % compute average time step size
    dt = dt/N;
end

% ------------------------------ END OF CODE ------------------------------
