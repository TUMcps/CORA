function traj = createTestSuite(sys, params, n_k, n_m, n_s, options)
% createTestSuite - Creates a trajectory array consisting of random
%   trajectories for the given system dynamics and parameters
%
% Syntax:
%    traj = createTestSuite(sys, params, n_k, n_m, n_s, options)
%    traj = createTestSuite(sys, params, n_k, n_m, n_s)
%
% Inputs:
%    sys - dynamical system
%    params - parameter defining the conformance problem
%       .R0 - initial state set
%       .U - input set
%    n_k - length per trajectory
%    n_m - number of trajectories
%    n_s - numbers of random simulations per trajectory
%    options - [optional] options specifying the sampling methods
%       .p_extr - fraction of simulations using extreme points
%           if not specified: 0
%       .inputCurve - array with form of the nominal input trajectories 
%                     ("rand", "randn", "bezier", "sinWave", "sigmoid") 
%                     for each input dimension
%           if not specified: "randn"
%       .inputSet - set of admissible parameters for generating nominal 
%                   input 
%           if not specified: sampling from default parameter sets
%       .stateSet - set of admissible nominal initial states
%           if not specified: nominal initial states are zero
%           
%
% Outputs:
%    traj - array of N_m trajectory objects
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Laura Luetzow, Zeqi Li
% Written:       27-August-2025            
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% settings
if nargin < 6
    options = {}; % use default settings
end

% probability of extreme samples
if ~isfield(options, 'p_extr')
    options.p_extr = 0;
end

% simulationi time
params.tFinal = n_k *sys.dt - sys.dt;

traj(n_m,1) = trajectory();
for m = 1:n_m
    % sample a random nominal input trajectory of the specified curve type
    u_m = aux_createCurve(sys.nrOfInputs,n_k,sys.dt,options);

    % sample a random nominal initial state vector
    if isfield(options, 'stateSet')
        x0_m = randPoint(options.stateSet, 1);
    else
        x0_m = zeros(dim(params.R0),1);
    end

    y_m = zeros(sys.nrOfOutputs,n_k,n_s);
    % sample the true input trajectory and initial state vector
    for s = 1:n_s
        if rand(1) < options.p_extr
            params.x0 = x0_m + randPoint(params.R0,1,'extreme');
        else
            params.x0 = x0_m + randPoint(params.R0,1);
        end
        if rand(1) < options.p_extr
            params.u = u_m + randPoint(params.U, n_k,'extreme');
        else
            params.u = u_m + randPoint(params.U, n_k);
        end

        % simulate the system dynamics to obtain the output trajectory
        [~,~,~,y_m(:,:,s)] = simulate(sys,params);
    end
    % create trajectory object for each nominal trajectory
    traj(m) = trajectory(u_m, x0_m, y_m, [], sys.dt, [], [], class(sys));
end
end


% Auxiliary functions -----------------------------------------------------

function u = aux_createCurve(n_u,n_k,dt,options)
% createCurve - generate curve from given parameters

tFinal = n_k*dt;
u = zeros(n_u,n_k);

% define curve type
if ~isfield(options, 'inputCurve')
    curveType = repmat("randn",n_u,1);
else
    curveType = options.inputCurve;
    if ~isstring(curveType)
        curveType = string(curveType);
    end
    if length(curveType) < n_u
        curveType = repmat(curveType,n_u,1);
    end
end

% loop through all input dimensions
for i = 1: n_u
    if isfield(options, 'inputSet') && iscell(options.inputSet)
        p_i = randPoint(inputSet{i}, 1);
    elseif isfield(options, 'inputSet')
        p_i = randPoint(inputSet, 1);
    else
        p_i = randPoint(interval([-10; 1; 1; 0; 0; 0],...
            [-1; 10; 1; 0; 2; 2]));
    end
    % minimum and maximum input with u_max>=u_min+1e-3
    u_min = p_i(1); 
    u_max = max(u_min+1e-3, p_i(2)); 
    % number of nonzero inputs (must be between 3 and n_k)
    n_nonzero = min(max(round(p_i(3)*n_k), 3), n_k);
    % start time step of nonzero inputs (must be smaller equal than n_k-n_nonzero+1)   
    k_start_max = n_k-n_nonzero+1;
    k_start = min(max(round(p_i(4)*k_start_max), 1), k_start_max); 
    % scaling parameters
    A1 = p_i(5); 
    A2 = p_i(6); 

    % generate time vector
    t = 1:n_nonzero;

    switch curveType(i)
        % generate input signals
        case "rand"
            % uniformly distributed
            u_i = u_min + (u_max-u_min) * rand(1,n_nonzero);

        case "randn"
            % Gaussian distributed
            u_i = A1 * randn(1,n_nonzero);

        case "bezier"
            % Bezier curve
            % control points
            P = [0 0; 0.25 A1; 0.75 A2; 1 0];  
            n = size(P,1)-1; % Degree of curve
            x = linspace(0,1,100);

            % generate input signals
            B = zeros(length(x),2);
            for j = 0:n
                % Bernstein polynomial
                B = B + (nchoosek(n,j) * (1-x).^(n-j) .* x.^j)' * P(j+1,:);
            end
            u_i = interp1(B(:,1), B(:,2), (t-1)/(n_nonzero-1), 'linear');

        case "sigmoid"
            % Signmoid curve
            % generate input signals;
            t_mid = t(round(n_nonzero/2));
            u_i = A1* 1./(1+A2/n_nonzero*exp(-t+t_mid));
            u_i = u_i - u_i(1); % start with zero

        case "sinWave"            
            % Sinus curve
            % generate input signals
            u_i = A1*sin(A2*(2*pi/n_nonzero)*(t-1));
    end
    u(i, k_start:k_start+n_nonzero-1) = u_i; 
    u(i, u(i,:) < u_min) = u_min;
    u(i, u(i,:) > u_max) = u_max;
end
end

% ------------------------------ END OF CODE ------------------------------
