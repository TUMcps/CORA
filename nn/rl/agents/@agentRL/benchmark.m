function [reachSet,reward,rewardAdvNaive,rewardAdvGrad] = benchmark(obj,env,perturbation,varargin)
% benchmark - benchmark the trained agentRL
%
% Syntax:
%   [reachSet,volumeLoss] = benchmark(obj,state,pertubation,initalOps)
%
% Inputs:
%   env - control environment
%   perturbation - current perturbation radius (scalar)
%   initialOps - inintial observation option for evaluation
% 
% Outputs:
%   reachSet - reachSet of the simulation 
%   volumeLossAction - volume loss of the action
%   columeLossState - volume loss of the system state
%
%                                                            
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: critic

% Authors:       Manuel Wendl
% Written:       04-November-2023
% Last update:   18-August-2024 (new options structure)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parse function arguments
narginchk(3,4)
initalOps = setDefaultValues('inf',varargin);
% Validate function arguments
inputArgsCheck({ ...
        {env, 'att', 'ctrlEnvironment'}; ...
        {perturbation, 'att', 'numeric'}; ...
        {initalOps, 'str', {'None','inf','sup'}}; ...
        })

obj.options.rl.noise = perturbation;   

if numel(obj.options.rl.noise) == 1
    % If pertubation radius is given
    noise = zonotope(zeros(obj.actor.nn.neurons_in,1),eye(obj.actor.nn.neurons_in)*obj.options.rl.noise);
else
    % custom pertubation generators are checked
    if size(obj.options.rl.noise,1) ~= obj.actor.nn.neurons_in
        throw(CORAerror("CORA:dimensionMismatch",obj.options.rl.noise,obj.actor.nn.neurons_in));
    else
        noise = zonotope(zeros(obj.actor.nn.neurons_in,1),obj.options.rl.noise);
        obj.options.rl.noise = mean(sum(abs(obj.options.rl.noise),2));
    end
end

% Reset Environment for set-based evaluation
env.options.rl.env.evalMode = 'set';
env.options.rl.env.initialOps = initalOps;
env.options.env.collisioCheckBool = false;

[env, observation] = env.reset();

% Instantiate output
reward = 0;
rewardAdvNaive = 0;
rewardAdvGrad = 0;
reachSet = [];

% Set based verification of adverserial attacks
isDone = false;
while ~isDone
    action = obj.actor.nn.evaluate_(observation+noise,obj.options.rl.actor,obj.actor.idxLayer);
    [env, nextObs, stepReward, isDone , visual] = env.step(action);  
    observation = nextObs;
    % Accumulate reward
    reward = reward + stepReward;
    % Store reachable set
    if isempty(reachSet)
        reachSet = visual;
    else
        reachSet =  add(reachSet,visual);   
    end
end

% Convert to point-based evaluation
env.options.rl.env.evalMode = 'point';
env.options.rl.env.initialOps = initalOps;
[env,observation] = env.reset();

% Store original actor train method in temp variable
original_attack = obj.options.rl.actor.nn.train.method;

% Adverserial Attack Naive
obj.options.rl.actor.nn.train.method = 'naive';

isDone = false;
while ~isDone
    advObs = obj.computeAdversarialAttack(observation);
    action = obj.actor.nn.evaluate_(advObs,obj.options.rl.actor,obj.actor.idxLayer);
    [env, nextObs, stepReward, isDone, ~] = env.step(action);  
    observation = nextObs;
    % Accumulate reward
    rewardAdvNaive = rewardAdvNaive + stepReward;
end


% Convert to point-based evaluation
env.options.rl.env.evalMode = 'point';
env.options.rl.env.initialOps = initalOps;
[env,observation] = env.reset();

% Adverserial Attack Naive
obj.options.rl.actor.nn.train.method = 'grad';

isDone = false;
while ~isDone
    advObs = obj.computeAdversarialAttack(observation);
    action = obj.actor.nn.evaluate_(advObs,obj.options.rl.actor,obj.actor.idxLayer);
    [env, nextObs, stepReward, isDone, ~] = env.step(action);  
    observation = nextObs;
    % Accumulate reward
    rewardAdvGrad = rewardAdvGrad + stepReward;
end

% Restore original train method
obj.options.rl.actor.nn.train.method = original_attack;
end 


% ------------------------------ END OF CODE ------------------------------
