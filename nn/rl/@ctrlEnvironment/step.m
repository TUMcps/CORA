function [obj, observation, reward, isDone, visual] = step(obj,action)
% step - one step of the ctrlEnvironment
%
% Syntax:
%   [observation, reward, isDone, visual] = step(obj,action)
%
% Inputs:
%   action - current action of the controller
% 
% Outputs:
%   env - updated environment
%   observation - end state of step
%   reward - reward for state and action
%   isDone - terminal flag
%   visual - visualisation data for rendering the environnment step
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ctrlEnvironment

% Authors:       Manuel Wendl
% Written:       03-November-2023
% Last update:   18-August-2024 (simplify step functionality by extracting from RLagent)
%                19-September-2024 (TL, rewrote reach call w/ params & options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({ ...
    {action, 'att', {'numeric','zonotope','gpuArray'}}; ...
    })

activeGPU = isgpuarray(action);

% case for Training Mode
if strcmp(obj.options.rl.env.evalMode,'point')
    % update time interval
    obj = obj.updateTimeInterval();
    
    if isa(action,'zonotope')
        throw(CORAerror("CORA:notSupported",'For point-based environement evaluation, no zononotope actions supported.'));
    end

    % Load simaction from GPU to CPU if necessary
    if activeGPU
        action = gather(action);
    end
    
    % simulate control time step
    ops.u = action;
    tVec = single(obj.options.rl.env.reach.tVec);
    
    if strcmp(obj.options.rl.env.solver,'ODE45')
        [t, x] = ode45(getfcn(obj.ctrlDynamics, ops), tVec, [obj.state;action]);
    elseif strcmp(obj.options.rl.env.solver,'Euler')
        t = obj.stepNum*obj.options.rl.env.dt;
        f = getfcn(obj.ctrlDynamics, ops);
        x = ([obj.state;action] + obj.options.rl.env.dt*f(t,[obj.state;action]))';
    else
        throw(CORAerror("CORA:notSupported",['No other solver ' ...
            'supported in the ctrlEnvironment than ODE45 and Euler.']));
    end

    % get reward function output
    if nargin(obj.rewardFun) == 1 % Stabilization tasks
        reward = obj.rewardFun(x);
    elseif nargin(obj.rewardFun) == 2 % Trajectory tasks
        reward = obj.rewardFun(x,t);
    else
        throw(CORAerror("CORA:specialError",['Reward Fun with more ' ...
            'than 2 input arguments is not supported!']))
    end

    % get next state from simulation
    observation = x(end, 1:(obj.ctrlDynamics.nrOfOutputs-obj.ctrlDynamics.nrOfInputs))';
    visual.t = t;
    visual.x = x;
    
    % Check if env is done 
    if obj.stepNum == obj.options.rl.env.maxSteps || (obj.collisionCheck(x(:,1:(obj.ctrlDynamics.nrOfOutputs-obj.ctrlDynamics.nrOfInputs)))&&obj.options.rl.env.collisionCheckBool)
        isDone = true;
    else
        isDone = false;
    end

else % case for set-based evaluation 
    % update time interval
    obj = obj.updateTimeInterval();
    
    if activeGPU
        action = zonotope(double(gather(action.c)),double(gather(action.G)));
    else
        % Cast action to double due to reach compatibility with sparse
        % double Hessian and third order tensors. 
        action = zonotope(double(action.c),double(action.G));
    end

    % append control input to initial set and cast state to double.
    obj.options.rl.env.reach.R0 = aux_appendAction(zonotope(double(obj.state.c),double(obj.state.G)), action);
    
    % get reachable set
    [R, res] = aux_reachability(obj.ctrlDynamics, obj.options.rl.env.reach);

    % get function output
    visual = R;
    reward = obj.rewardFun(R);
    if obj.stepNum < obj.options.rl.env.maxSteps && ~res
        observation = R.timePoint.set{1};
        reward = -1e8;
        isDone = true;
        return;
    else
        R = project(R, 1:(obj.ctrlDynamics.nrOfOutputs-obj.ctrlDynamics.nrOfInputs));
    end
    
    % Get observation from last state of reachability analysis and cast
    % back to single.
    observation = zonotope(single(R.timePoint.set{end}.c),single(R.timePoint.set{end}.G));

    if obj.stepNum == obj.options.rl.env.maxSteps || (obj.collisionCheck(R)&&obj.options.rl.env.collisionCheckBool)
        isDone = true;
    else
        isDone = false;
    end
end

% update environment
obj.state = observation;
obj.stepNum = obj.stepNum + 1;
end


% Auxiliary functions -----------------------------------------------------

function [R, res] = aux_reachability(sys, options)
% compute the reachable set of the controlled system (without performing
% option checks)

% rewrite options to params (already validated)
params = struct;
fields = {'tStart','tFinal','R0','U','uTrans','factor','tVec'};
for i = 1:numel(fields)
    field = fields{i};
    if isfield(options,field)
        params.(field) = options.(field);
        options = rmfield(options,field);
    end
end

% call reach
[R,res] = reach(sys,params,options);

end

% append initial zonotope and action
function R0 = aux_appendAction(X,U)
G0 = [X.G, zeros(size(X.c, 1), size(U.G, 2)-size(X.G, 2)); U.G];
R0 = zonotope(cat(1,X.c,U.c),G0);
end

% ------------------------------ END OF CODE ------------------------------
