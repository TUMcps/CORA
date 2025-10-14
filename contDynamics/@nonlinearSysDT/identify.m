function sys = identify(varargin)
% identify - Identify a nonlinear discrete-time system from trajectory data
%
% Syntax:
%    sys = nonlinearSysDT.identify(traj)
%    sys = nonlinearSysDT.identify(traj,options)
%    sys = nonlinearSysDT.identify(x,t)
%    sys = nonlinearSysDT.identify(x,t,options)
%    sys = nonlinearSysDT.identify(x,t,u)
%    sys = nonlinearSysDT.identify(x,t,u,options)
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
%       -.dt:       sampling time for the discrete-time system. The default
%                   value is the average time step size from the data. 
%       -.alg:      algorithm that is used ('dmd': Dynamic Mode 
%                   Decomposition [1], 'opt': optimization, 'sparse': SINDy 
%                   approach [2], 'nested': nested template functions).
%                   The default value is 'dmd'.
%       -.basis:    template/basis functions for the dynamics 
%                   ('quad': quadratic terms, 'cubic': cubic terms, 
%                    'trig': trigonometric functions, 'all': all together).
%                   The default value is 'trig'. (all algorithms except
%                   'nested' and only if options.phi is not specified)
%       -.phi:      user-defined template functions provided as a function
%                   handle, e.g. options.phi = @(x,u) [x(1);x(2)*u(1)] 
%                   (all algorithms except 'nested')
%       -.inpAff:   boolean specifying if the identified model should be 
%                   input-affine (options.inpAff = true) or not. The 
%                   default value is options.inpAff = true. 
%
% Outputs:
%    sys - identified nonlinearSysDT object
%
% Example: 
%    f = @(x,u) [0.99*x(1) + 0.2*x(2); ...
%                -0.1*x(1) + 0.5*x(2)/(1+x(2)^2)];
%    dt = 1;
%    sysOrig = nonlinearSysDT(f,dt);
%
%    simOpts.x0 = [-8; 7];
%    simOpts.tFinal = 10;
%    [t,x] = simulate(sysOrig,simOpts);
%
%    sys = nonlinearSysDT.identify(x,t);
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
%    [2] S.L. Brunton and et al. "Discovering governing equations from data 
%         by sparse identification of nonlinear dynamical systems, 2016
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys/identify

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
            phi = aux_nonlinearBasisFunctions(traj,options);
            sys = priv_identifyDMD(traj,phi);

        case 'opt'
            phi = aux_nonlinearBasisFunctions(traj,options);
            sys = priv_identifyOpt(traj,phi);

        case 'sparse'
            phi = aux_nonlinearBasisFunctions(traj,options);
            sys = priv_identifySparse(traj,phi);

        case 'nested'
            sys = priv_identifyNested(traj);
    end
end


% Auxiliary functions -----------------------------------------------------

function options = aux_checkOptions(options,traj)
% check the algorithm settings provided by the user

    % check if there are any redundant options specified
    if ~isempty(options)
        if isfield(options,'alg') && strcmp(options.alg,'nested')
            redundantOptions(options,{'alg','dt'});
        elseif isfield(options,'phi')
            redundantOptions(options,{'alg','dt','phi'});
        else
            redundantOptions(options,{'alg','dt','phi','basis','inpAff'});
        end
    end

    % default values
    if ~isfield(options,'alg')
        options.alg = 'dmd';
    end

    if ~isfield(options,'basis')
        options.basis = 'trig';
    end

    if ~isfield(options,'inpAff')
        options.inpAff = true;
    end

    if ~isfield(options,'dt')
        options.dt = aux_averageTimeStepSize(traj);
    end

    % check user input (sampling time)
    if options.dt <= 0 || ~isscalar(options.dt)
        throw(CORAerror('CORA:wrongFieldValue','options.dt','double > 0'));
    end

    % check user input (specified algorithm)
    if ~ismember(options.alg,{'dmd','opt','sparse','nested'})
        throw(CORAerror('CORA:wrongFieldValue','options.alg', ...
                                        {'dmd','opt','sparse','nested'}));
    end

    % check user input (basis functions)
    if ~ismember(options.basis,{'quad','cubic','trig','all'})
        throw(CORAerror('CORA:wrongFieldValue','options.basis', ...
                                        {'quad','cubic','trig','all'}));
    end

    % check user input (custom template function)
    if isfield(options,'phi')
        x = sym('x',[size(traj(1).x,1),1]);
        u = sym('u',[min(1,size(traj(1),1)),1]); %TO-DO: traj(1).u???
        try
            tmp = options.phi(x,u);
            if size(tmp,2) ~= 1
                throw(CORAerror('CORA:specialError', ...
                      'options.phi: function handle has to return a vector'));
            end
        catch
            throw(CORAerror('CORA:specialError', ...
                      'options.phi: error in function handle evaluation'));
        end
    end

    % check user input (input affine)
    if ~islogical(options.inpAff)
        throw(CORAerror('CORA:wrongFieldValue','options.inpAff','boolean'));
    end
end

function phi = aux_nonlinearBasisFunctions(traj,options)
% construct nonlinear template functions for the nonlinear system

    % check if custom user-defined template functions are provided
    if isfield(options,'phi')
        phi = options.phi; return;
    end

    % construct basic functions
    n = size(traj(1).x,1);
    x = sym('x',[n,1]);

    temp = num2cell(x);
    poly = temp;

    if ismember(options.basis,{'trig','all'})
        poly = [poly; cellfun(@(x) cos(x),temp,'UniformOutput',false)];
        poly = [poly; cellfun(@(x) sin(x),temp,'UniformOutput',false)];
    end

    % add inputs
    if ~isempty(traj(1).u)
        m = size(traj(1).u,1);
        u = sym('u',[m,1]);
        if options.inpAff
            poly = [poly; num2cell(u)];
        end
    else
        u = sym('u',[1,1]);
    end

    % add all quadratic combinations of the basis functions
    if ismember(options.basis,{'quad','trig'})

        poly_ = poly;
    
        for i = 1:length(poly_)
            for j = i:length(poly_)
                poly{end+1} = poly_{i} * poly_{j};
            end
        end
    end

    % add all cubic combinations of the basis functions 
    if ismember(options.basis,{'cubic','all'})
        
        poly_ = poly;
    
        for i = 1:length(poly_)
            for j = i:length(poly_)
                for k = j:length(poly_)
                    poly{end+1} = poly_{i} * poly_{j} * poly_{k};
                end
            end
        end
    end

    % add constant input
    poly = cellfun(@(x) x,poly);
    poly = [poly;sym(1)];

    % add input
    if options.inpAff && ~isempty(traj(1).u)
        poly = [poly;u];
    end

    % convert to function handle
    phi = matlabFunction(poly,'Vars',{x,u});
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
