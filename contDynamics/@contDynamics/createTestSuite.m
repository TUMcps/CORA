function testSuite = createTestSuite(sys, params, n_k, n_m, n_s, options)
% createTestSuite - Creates a testSuite consisting of test cases for
%   the given system dynamics
%
% Syntax:
%    testSuite = createTestSuite(sys, params, n_k, n_m, n_s, options)
%    testSuite = createTestSuite(sys, params, n_k, n_m, n_s)
%
% Inputs:
%    sys - dynamical system
%    params - parameter defining the conformance problem
%       .R0 - initial state set
%       .U - input set
%    n_k - length per test case
%    n_m - number of test cases
%    n_s - numbers of random simulations per test case
%    options - [optional] options specifying the sampling methods
%       .p_extr - fraction of simulations using extreme points
%           if not specified: 0
%       .p_biased - fraction of simulation from biased uncertainty set
%           if not specified: 0
%       .inputCurve - array with form of the nominal input trajectories 
%                     ("rand", "randn", "bezier", "sinWave", "sigmoid") 
%                     for each input dimension
%           if not specified: "randn"
%       .inputParameters - n_p x n_u array with the curve parameters
%           if not specified: sampling from the inputSet curve 
%       .inputSet - set of admissible parameters for nominal input 
%                   (or as cell array with for each input dimension)
%           if not specified: sampling from default parameter sets
%       .inputFactor - vector with a multiplicative scaling factor
%                      for each input dimension
%           if not specified: no scaling
%       .stateSet - set of admissible nominal initial states
%           if not specified: nominal initial states are zero
%       .contInput - boolean specifying if continuous input is required
%           if not specified: true
%           
%
% Outputs:
%    testSuite - cell array of N_m testCase objects
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Laura Luetzow, Zeqi Li
% Written:       28-February-2024             
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 6
    options = {}; % use default settings
end
if ~isfield(options, 'p_extr')
    options.p_extr = 0;
end
if ~isfield(options, 'p_biased')
    options.p_biased = 0;
elseif ~isa(params.R0, 'zonotope') || ~isa(params.U, 'zonotope')
    throw(CORAerror('CORA:notSupported',"Creation of biased test cases is only implemented for zonotopic uncertainty sets."))
else
    U_biased = zonotope(center(params.U)+ sum(0.5*generators(params.U),2), 0.5*generators(params.U));
end

% figure; hold on; 
params.tFinal = n_k *sys.dt - sys.dt;

if isa(sys,'linearSysDT') || isa(sys,'linearARX')
    u=[];
    x0 = [];
    y = [];
else
    testSuite = cell(n_m,1);
end
for m = 1:n_m
    % sample a random nominal input trajectory of the specified curve type
    u_m = aux_createCurve(sys.nrOfInputs,n_k,sys.dt,options);

    % sample a random nominal initial state vector
    if isfield(options, 'stateSet')
        x0_m = randPoint(options.stateSet, 1);
    else
        x0_m = zeros(dim(params.R0),1);
    end

    y_m = zeros(n_k,sys.nrOfOutputs,n_s);
    % sample the true input trajectory and initial state vector
    for s = 1:n_s
        if rand(1) < options.p_biased
            U = U_biased;
        else
            U = params.U;
        end
        R0 = params.R0;
        if rand(1) < options.p_extr
            params.x0 = x0_m + randPoint(R0,1,'extreme');
        else
            params.x0 = x0_m + randPoint(R0,1);
        end
        if rand(1) < options.p_extr
            params.u = u_m + randPoint(U, n_k,'extreme');
        else
            params.u = u_m + randPoint(U, n_k);
        end

        % simulate the system dynamics to obtain the output trajectory
        [~,~,~,y_m(:,:,s)] = simulate(sys,params);
    end
    if n_m > 1 && (isa(sys,'linearSysDT') || isa(sys,'linearARX'))
        % create one test case by concatenating all trajectories in the 
        % third dimension
        x0 = cat(3,x0,repmat(x0_m,1,1,n_s));
        u = cat(3,u,repmat(u_m',1,1,n_s));
        y = cat(3,y,y_m);
    else
        % create testCase for each nominal trajectory
        testSuite{m} = testCase(y_m, u_m', x0_m, sys.dt, class(sys));
    end
end

if n_m > 1 && (isa(sys,'linearSysDT') || isa(sys,'linearARX'))
    testSuite{1} = testCase(y, u, x0, sys.dt, class(sys));
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

% define curve parameters or parameter set for sampling
default_inputSet = false;
if isfield(options, 'inputParameters') 
    % curve parameters are given
    p = options.inputParameters;

else 
    % curve parameters have to be sampled
    if ~isfield(options, 'inputSet')
        % sampling from the default parameter set
        default_inputSet = true; 
        d_min = 0; % minimum time delay
        d_max = 1/4*tFinal; % maximum time delay

    elseif iscell(options.inputSet)
        % sampling from different given parameter sets
        inputSet = options.inputSet; 

    else
        % sampling from the same given parameter set
        inputSet = repmat({options.inputSet},n_u,1);
    end
end

for i = 1: n_u
    if isfield(options, 'inputParameters')
        p_i = p(:,i);
    elseif isfield(options, 'inputSet')        
        p_i = randPoint(inputSet{i}, 1);
    end
    switch curveType(i)
        case "rand"
            if default_inputSet
                p_i = randPoint(interval([d_min;-1;3/4*tFinal],[d_max;1;tFinal]),1);
            end

            % load parameters of test case-------------------------------------
            delay = max(p_i(1),0); % timeDelay
            A = p_i(2); % amplitude
            tZero = max(min(p_i(3),tFinal), delay+dt); % time after which input is equal to zero

            % generate input signals-------------------------------------------
            t = max(ceil( delay /dt)*dt,dt) : dt : min( floor( (tZero+delay) /dt)*dt, tFinal);
            u(i,max(ceil( delay /dt),1): min(floor( (tZero+delay) /dt),n_k)) = ...
                A * rand(length(t),1);

        case "randn"
            if default_inputSet
                p_i = randPoint(interval([d_min;-1;3/4*tFinal],[d_max;1;tFinal]),1);
            end

            % load parameters of test case-------------------------------------
            delay = max(p_i(1),0); % timeDelay
            A = p_i(2); % amplitude
            tZero = max(min(p_i(3),tFinal), delay+dt); % time after which input is equal to zero

            % generate input signals-------------------------------------------
            t = max(ceil( delay /dt)*dt,dt) : dt : min( floor( (tZero+delay) /dt)*dt, tFinal);
            u(i,max(ceil( delay /dt),1): min(floor( (tZero+delay) /dt),n_k)) = ...
                A * randn(length(t),1);

        case "bezier"
            if default_inputSet
                p_i = randPoint(interval([d_min;-1;3/4*tFinal;0;0],[d_max;1;tFinal;8;8]),1);
            end

            % load parameters of test case-------------------------------------
            delay = max(p_i(1),0); % timeDelay
            A = p_i(2); % slipCompensation
            tZero = max(min(p_i(3),tFinal), delay+dt);
            p1 = p_i(4); % \in [0,8]
            p2 = p_i(5); % \in [0,8]

            % generate input signals-------------------------------------------
            t = max(ceil( delay /dt)*dt,dt) : dt : min( floor( (tZero+delay) /dt)*dt, tFinal);
            u(i,max(ceil( delay /dt),1): min(floor( (tZero+delay) /dt),n_k)) = ...
                A*( 30/tZero + (6*p1-120)/tZero^2*(t-delay) + ...
                (90-9*p1+3*p2)/tZero^3*(t-delay).^2 );

        case "sigmoid"
            if default_inputSet
                p_i = randPoint(interval([d_min;-1;-1; -1; 1],[d_max;1;1;1;10]),1);
            end

            % load parameters of test case-------------------------------------
            delay = max(p_i(1),0); % timeDelay
            A = p_i(2); % slipCompensation
            A1 = p_i(3); % \in [2,4]
            A2 = p_i(4); % \in [-2,4]
            p3 = p_i(5); % \in [4,12]

            % generate input signals-------------------------------------------
            t = max(ceil( delay /dt)*dt,dt) : dt : tFinal;
            u(i,max(ceil( delay /dt),1):n_k) = ...
                A*(A1*p3.*exp(-p3.*((t-delay)-0.5))./(1+exp(...
                -p3.*((t-delay)-0.5))).^2 + A2*p3.*exp(-p3*((t-delay)-2))./...
                (1+exp(-p3.*((t-delay)-2))).^2);

        case "sinWave"
            if default_inputSet
                p_i = randPoint(interval([d_min;-1;tFinal/10; 0],[d_max;1;10*tFinal;n_k]),1);
            end

            % load parameters of test case-------------------------------------
            delay = max(p_i(1),0); % timeDelay
            A = p_i(2); % amplitude
            T = 10*max(dt,p_i(3))*dt; 
            Ts = max(0,p_i(4))*dt; % duration for keeping a constant input

            % generate input signals-------------------------------------------
            t1 = max(ceil(delay/dt)*dt,dt) : dt : min(floor((T/4+delay)/dt),n_k)*dt;
            t2 = ceil((T/4+Ts+delay)/dt)*dt : dt : min(floor((T/2+Ts+delay)/dt),n_k)*dt ;

            u(i,max(ceil(delay/dt),1):min(floor((T/4+delay)/dt),n_k)) = ...
                A*cos((2*pi/T)*(t1-delay));
            u(i,ceil((T/4+Ts+delay)/dt):min(floor((T/2+Ts+delay)/dt),n_k)) = ...
                A*cos((2*pi/T)*(t2-Ts-delay));
    end

    if isfield(options, 'inputFactor')
        u(:,i) = options.inputFactor(i) * u(:,i);
    end
end
if ~isfield(options, 'contInput') || options.contInput
    u = cumsum(u,2);
end
end

% ------------------------------ END OF CODE ------------------------------
