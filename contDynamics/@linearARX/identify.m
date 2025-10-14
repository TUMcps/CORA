function sys = identify(varargin)
% identify - Identifies a linear ARX system from trajectory data
%
% Syntax:
%    sys = linearARX.identify(traj)
%    sys = linearARX.identify(traj,options)
%    sys = linearARX.identify(y,t)
%    sys = linearARX.identify(y,t,options)
%    sys = linearARX.identify(y,t,u)
%    sys = linearARX.identify(y,t,u,options)
%
% Inputs:
%    traj - object of class "trajectory" storing the trajectory data
%    y - cell-array storing the output of the simulated trajectories, where
%        each trajectory is a matrix of dimension [n,N]
%    t - cell-array storing the time points for the simulated trajectories
%        as vectors of dimension [1,N]
%    u - cell-array storing the inputs for the simulated trajectories
%        as vectors of dimension [m,N-1]
%    options - algorithm options for system identification
%
%       -.dt:   sampling time for the discrete-time system. The default
%               value is the average time step size from the data.
%       -.p:    number of past time steps. The default value is 2.
%       -.alg:  algorithm that is used ('dmd': Dynamic Mode 
%               Decomposition [1], 'opt': optimization). The default value
%               is 'dmd'.
%
% Outputs:
%    sys - identified linearARX object
%
% Example: 
%    dt = 0.1;
%    A_bar = {[-0.4 0.6; 0.6 -0.4];[0.1 0; 0.2 -0.5]};
%    B_bar = {[0; 0];[0.3; -0.7];[0.1; 0]};
%    linARX = linearARX(A_bar,B_bar,dt);
%
%    params.y0 = [[0;0] [0.1;0.2]];
%    params.tFinal = 4;
%    params.u = [repmat([0.1 0.05 0.05 -0.1],[1,10]),0.1];
%
%    [t,~,~,y] = simulate(linARX,params);
%
%    sys = linearARX.identify(y,t,params.u);
%
%    [t_,~,~,y_] = simulate(sys,params);
%
%    figure; hold on; box on;
%    plot(t,y(1,:));
%    plot(t_,y_(1,:),'LineStyle','--');
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
% Written:       30-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % replace the states of trajectory with the outputs if outputs exist
    varargin = aux_output2state(varargin);

    % parse the input arguments and bring trajectory data to suitable form
    [traj,options] = checkDataIdentification(varargin{:});

    % check the algorithm options
    options = aux_checkOptions(options,traj);

    % convert data to uniform time step size
    traj = uniformTimeStepSize(traj,options.dt);

    % call the selected algorithm
    switch options.alg

        case 'dmd'
            sys = priv_identifyDMD(traj,options.p);

        case 'opt'
            sys = priv_identifyOpt(traj,options.p);
    end
end


% Auxiliary functions -----------------------------------------------------

function options = aux_checkOptions(options,data)
% check the algorithm settings provided by the user

    % default values
    if ~isfield(options,'alg')
        options.alg = 'dmd';
    end

    if ~isfield(options,'p')
        options.p = 2;
    end

    if ~isfield(options,'dt')
        options.dt = aux_averageTimeStepSize(data);
    end

    % check user input
    if options.dt <= 0 || ~isscalar(options.dt)
        throw(CORAerror('CORA:wrongFieldValue','options.dt','double > 0'));
    end

    if options.p <= 0 || ~isscalar(options.p) || mod(options.p,1) ~= 0
        throw(CORAerror('CORA:wrongFieldValue','options.p','integer > 0'));
    end

    if ~ismember(options.alg,{'dmd','opt'})
        throw(CORAerror('CORA:wrongFieldValue','options.alg', ...
                                                        {'dmd','opt'}));
    end

    % check if there are any redundant options specified
    redundantOptions(options,{'alg','dt','p'});
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

function inpArgs = aux_output2state(inpArgs)
% replace the states of the trajectory with the outputs if outputs exist

    if isa(inpArgs{1},'trajectory')

        % replace states with outputs for all trajectories
        traj = inpArgs{1}; flag = zeros(length(traj),1);

        for i = 1:length(traj)
            if ~isempty(traj(i).y)
                traj(i) = trajectory(traj(i).u,traj(i).y,[],traj(i).t);
                flag(i) = 1;
            end
        end

        inpArgs{1} = traj;

        % check if the provided data format is correct
        if length(unique(flag)) > 1
            throw(CORAerror('CORA:wrongValue','first', ...
                  'trajectory object with consistent states and outputs'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
