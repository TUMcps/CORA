function sys = identify(varargin)
% identify - Identifies a linear system from trajectory data
%
% Syntax:
%    sys = linearSys.identify(traj)
%    sys = linearSys.identify(traj,options)
%    sys = linearSys.identify(x,t)
%    sys = linearSys.identify(x,t,options)
%    sys = linearSys.identify(x,t,u)
%    sys = linearSys.identify(x,t,u,options)
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
%       -.alg:  algorithm that is used ('dmd': Dynamic Mode 
%               Decomposition [1], 'opt': optimization)
%
% Outputs:
%    sys - identified linear system object
%
% Example: 
%    A = [-0.7 -2; 2 -0.7];
%    sysOrig = linearSys(A);
%
%    simOpts.x0 = [10; 5];
%    simOpts.tFinal = 5;
%    [t,x] = simulate(sysOrig,simOpts);
%
%    sys = linearSys.identify(x,t);
%
%    [t_,x_] = simulate(sys,simOpts);
%
%    figure; hold on; box on;
%    plot(x(1,:),x(2,:));
%    plot(x_(1,:),x_(2,:));
%
% References:
%    [1] J.L. Proctor and et al. "Dynamic mode decomposition with control" 
%         SIAM Journal on Applied Dynamical Systems 15, 1 (2016), 142–161
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

    % parse the input arguments and bring trajectory data to suitable form
    [traj,options] = checkDataIdentification(varargin{:});

    % check the algorithm options
    options = aux_checkOptions(options);

    % call the selected algorithm
    switch options.alg

        case 'dmd'
            sys = priv_identifyDMD(traj);

        case 'opt'
            sys = priv_identifyOpt(traj);
    end
end


% Auxiliary functions -----------------------------------------------------

function options = aux_checkOptions(options)
% check the algorithm settings provided by the user

    % default values
    if ~isfield(options,'alg')
        options.alg = 'dmd';
    end

    % check user input
    if ~ismember(options.alg,{'dmd','opt'})
        throw(CORAerror('CORA:wrongFieldValue','options.alg', ...
                                                        {'dmd','opt'}));
    end

    % check if there are any redundant options specified
    redundantOptions(options,{'alg'});
end

% ------------------------------ END OF CODE ------------------------------
