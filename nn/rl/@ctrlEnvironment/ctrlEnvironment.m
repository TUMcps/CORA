classdef ctrlEnvironment 
% ctrlEnvironment - reinforcement learning control environment
%   The environment defines the transition function to the next observation
%
% Syntax:
%   obj = ctrlEnvironment(sysDynamics,rewardFun,collisionCheck,varargin)
%
% Inputs:
%   sysDynamics - nonlinearSys
%   rewardFun - reward function 
%   collisionCheck - collision check 
%   options.rl - reinforcement learning options 
%       .env - environment options
%           .x0: [-1,1] (default) Initial observation interval
%           .initialOps: 'uniform' (default) Sampling option for x0.
%           .evalMode: 'point' (default) Evaluation mode for ctrEnv.
%           .collisionCheckBool: true (default) Boolean for collision
%               check function evaluation.
%           .dt: .1 (default): Control time step between new control
%               signals.
%           .timeStep: .01 (default) Simulation timestep of environment.
%           .maxSteps: 30 (default) Maximum number of steps of env.
%           .solver: 'ODE45 (default) Solver for simualtion of transition
%               function.
%           .reach - options for reachable alg.
%               .alg: 'lin' (default) Reachability algorithm.
%               .zonotopeOrder: 300 (default) Maximum zonotope order.;
%
% Outputs:
%   obj - generated ctrlEnvironment
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ctrlEnvironment

% Authors:       Manuel Wendl
% Written:       22-October-2023
% Last update:   18-August-2024 (parse settings in constructor not in reset)
%                18-September-2024 (TL, updated validateOptions call)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    properties
        state
        stepNum
        ctrlDynamics
        rewardFun
        collisionCheck
        options
    end

    methods
        % constuctor
        function obj = ctrlEnvironment(sysDynamics,rewardFun,collisionCheck,varargin)
            % Parse function arguments
            narginchk(3,4)
            options = setDefaultValues(struct(),varargin);

            inputArgsCheck({ ...
                {sysDynamics, 'att','nonlinearSys'}, ...
                {rewardFun, 'att','function_handle'}, ...
                {collisionCheck, 'att','function_handle'}, ...
                {options, 'att','struct'}, ...
                })

            obj.options = aux_validateEnvOptions(options,sysDynamics);

            sysCL = aux_computeProperties(sysDynamics);
            obj.ctrlDynamics = sysCL; 

            obj.rewardFun = rewardFun;
            obj.collisionCheck = collisionCheck;

            obj.stepNum = 1;

            [params,ops] = obj.setDefaultReach(zonotope(obj.options.rl.env.x0),obj.options);
            [obj,~,obj.options.rl.env.reach] = obj.parseSettings(params,ops);
        end
    end

    methods (Access = private)

        % extract default reachability options
        function [params, ops] = setDefaultReach(obj,state,options)
            ops = options.rl.env.reach;
            ops.timeStep = options.rl.env.timeStep;

            params.tStart = (obj.stepNum-1)*options.rl.env.dt;
            params.tFinal = obj.stepNum*options.rl.env.dt;
            params.R0 = state;
        end

        % parse reachability settings
        function [obj, params, ops] = parseSettings(obj,params,ops)
            % check if the algorithm settings provided by the user are correct
            params.R0 = cartProd(params.R0, zeros(obj.ctrlDynamics.nrOfInputs, 1));
            [params,ops] = validateOptions(obj.ctrlDynamics, params, ops, 'FunctionName','reach');
            ops = params2options(params,ops);
            ops.R0 = project(ops.R0, 1:(obj.ctrlDynamics.nrOfOutputs-obj.ctrlDynamics.nrOfInputs));

            % obtain factors for initial state and input solution time step
            r = ops.timeStep;
            for i = 1:(ops.taylorTerms + 1)
                ops.factor(i) = r^(i) / factorial(i);
            end

            % check if splitting is turned off
            if ~all(isinf(ops.maxError))
                throw(CORAerror('CORA:notSupported',...
                    ['Splitting reachable sets is not supported for neural', ...
                    ' network controlled systems!']));
            end

            % pre-compute derivatives
            ops.verbose = false;
            derivatives(obj.ctrlDynamics, ops);
        end

        % update Time interval
        function obj = updateTimeInterval(obj)
            % update time interval
            obj.options.rl.env.reach.tStart = (obj.stepNum-1)*obj.options.rl.env.dt;
            obj.options.rl.env.reach.tFinal = obj.stepNum*obj.options.rl.env.dt;

            tVec = obj.options.rl.env.reach.tStart:obj.options.rl.env.timeStep:obj.options.rl.env.reach.tFinal;
            if tVec(end) ~= obj.options.rl.env.reach.tFinal
                % add tFinal if sampling time and time horizon don't match
                % resulting in a partial time step at the end.
                tVec(end+1) = obj.options.rl.env.reach.tFinal;
            end
            obj.options.rl.env.reach.tVec = tVec;
        end
    end 
end


% Auxiliary functions -----------------------------------------------------

% set default values for the DDPGagent
function options = aux_validateEnvOptions(options,sysDynamics)

persistent defaultEnvFields
if isempty(defaultEnvFields)
    defaultEnvFields = {
        'x0', interval(-ones(sysDynamics.nrOfStates,1),ones(sysDynamics.nrOfStates,1));
        'initialOps', 'uniform';
        'evalMode', 'point';
        'collisionCheckBool', true;
        'dt', .1;
        'timeStep', .01;
        'maxSteps', 30;
        'solver', 'ODE45';
        'reach', struct();
        };
end

persistent defaultEnvReachFields
if isempty(defaultEnvReachFields)
    defaultEnvReachFields = {
       'alg', 'lin';
       'tensorOrder', 3;
       'taylorTerms', 4;
       'zonotopeOrder', 200;
       'errorOrder', 10;
       'intermediateOrder', 50;
    };
end

if ~isfield(options,'rl')
    options.rl = struct;
end

% check if any env options are given
if ~isfield(options.rl,'env')
    options.rl.env = struct;
end

% set default env options if required
for i=1:size(defaultEnvFields, 1)
    field = defaultEnvFields{i, 1};
    if ~isfield(options.rl.env, field)
        fieldValue = defaultEnvFields{i, 2};
        if isa(fieldValue, "function_handle")
            fieldValue = fieldValue(options);
        end
        options.rl.env.(field) = fieldValue;
    end
end

% set default env reach options if required
for i=1:size(defaultEnvReachFields, 1)
    field = defaultEnvReachFields{i, 1};
    if ~isfield(options.rl.env.reach, field)
        fieldValue = defaultEnvReachFields{i, 2};
        if isa(fieldValue, "function_handle")
            fieldValue = fieldValue(options);
        end
        options.rl.env.reach.(field) = fieldValue;
    end
end

% Check env fields
if CHECKS_ENABLED
    structName = inputname(1);
    aux_checkFieldClass(options.rl.env,'x0',{'interval'},structName);
    aux_checkFieldStr(options.rl.env,'initialOps',{'uniform', 'symmetric', 'None', 'inf', 'sup', 'set'},structName);
    aux_checkFieldStr(options.rl.env,'evalMode',{'point', 'set'},structName);
    aux_checkFieldNumericDefInterval(options.rl.env,'dt',interval(0,inf),structName);
    aux_checkFieldNumericDefInterval(options.rl.env,'timeStep',interval(1e-8,options.rl.env.dt),structName);
    aux_checkFieldNumericDefInterval(options.rl.env,'maxSteps',interval(1,inf),structName);
    aux_checkFieldStr(options.rl.env,'solver',{'ODE45', 'Euler'},structName);
    aux_checkFieldStr(options.rl.env.reach,'alg',{'lin'},structName);
    aux_checkFieldNumericDefInterval(options.rl.env.reach,'zonotopeOrder',interval(size(options.rl.env.x0.inf,1),inf),structName);
end

end

function aux_checkFieldStr(optionsenv, field, admissibleValues, structName)
fieldValue = optionsenv.(field);
if ~(isa(fieldValue, 'string') || isa(fieldValue, 'char')) || ...
        ~ismember(fieldValue, admissibleValues)
    throw(CORAerror('CORA:wrongFieldValue', ...
        aux_getName(structName, field), admissibleValues))
end
end

function aux_checkFieldClass(optionsenv, field, admissibleClasses, structName)
if ~ismember(class(optionsenv.(field)), admissibleClasses)
    throw(CORAerror('CORA:wrongFieldValue', ...
        aux_getName(structName, field), admissibleClasses))
end
end

function aux_checkFieldNumericDefInterval(optionsenv, field, I, structName)
    if ~contains_(I,optionsenv.(field),'exact',eps)
        throw(CORAerror('CORA:outOfDomain', ...
            aux_getName(structName, field),"ValidDomain",I))
    end
end

function msg = aux_getName(structName, field)
    msg = sprintf("%s.nn.%s", structName, field);
end

function sys = aux_computeProperties(sys)
% compute properties of neurNetContrSys object

    n = sys.nrOfStates; m = sys.nrOfInputs;
    % instantiate closed-loop system
    f = @(x, u) [sys.mFile(x(1:n), [x(n+1:n+m); u]); zeros(m, 1)];
    name = [sys.name, 'Controlled'];
    sys = nonlinearSys(name, f, n+m, max(1, sys.nrOfInputs-m));


end

% ------------------------------ END OF CODE ------------------------------
