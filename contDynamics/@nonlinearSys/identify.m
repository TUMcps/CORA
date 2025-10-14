function sys = identify(varargin)
% identify - Identifies a nonlinear system from trajectory data
%
% Syntax:
%    sys = nonlinearSys.identify(traj)
%    sys = nonlinearSys.identify(traj,options)
%    sys = nonlinearSys.identify(x,t)
%    sys = nonlinearSys.identify(x,t,options)
%    sys = nonlinearSys.identify(x,t,u)
%    sys = nonlinearSys.identify(x,t,u,options)
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
%    sys - identified nonlinearSys object
%
% Example: 
%    sigma = 10; beta = 8/3; rho = 28;
%    f = @(x,u) [sigma*(x(2) - x(1)); ...
%                x(1)*(rho - x(3)) - x(2); ...
%                x(1)*x(2)- beta*x(3)];
%    lorenz = nonlinearSys(f);
%
%    simOpts.x0 = [-8; 7; 27];
%    simOpts.tFinal = 10;
%    [t,x] = simulate(lorenz,simOpts);
%
%    sys = nonlinearSys.identify(x,t);
%
%    [t_,x_] = simulate(sys,simOpts);
%
%    figure; hold on; box on;
%    plot3(x(1,:),x(2,:),x(3,:));
%    plot3(x_(1,:),x_(2,:),x_(3,:));
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
% See also: linearSys/identify

% Authors:       Niklas Kochdumper
% Written:       10-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse the input arguments and bring trajectory data to suitable form
    [traj,options] = checkDataIdentification(varargin{:});

    % check the algorithm options
    options = aux_checkOptions(options,traj);

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
            redundantOptions(options,{'alg'});
        elseif isfield(options,'phi')
            redundantOptions(options,{'alg','phi'});
        else
            redundantOptions(options,{'alg','phi','basis','inpAff'});
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
        u = sym('u',[min(1,size(traj(1),1)),1]); %TO-DO: traj(1).u??
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

% ------------------------------ END OF CODE ------------------------------
