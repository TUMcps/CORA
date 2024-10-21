classdef agentDDPG < agentRL
% agentDDPG - deep deterministic policy gradient agentRL
%
% Syntax:
%   obj = agentDDPG(actorNN,criticNN,varargin)
%
% Inputs:
%   actorNN - neural network
%   criticNN - neural network 
%   options.rl - options for RL:
%        .gamma: .99 (default) Discount factor applied to future rewards 
%               during training.       
%        .tau: .005 (default) Smoothing factor for target network updates
%        .expNoise: .2 (default) Exploration noise std. deviation.
%        .expNoiseType: 'OU' (default) Exploration noise type specified as 
%               an OrnsteinUhlenbeck (OU) or Gaussian.
%        .expDecayFactor: 1 (default) Exploration noise decay factor.
%               Linearly decreasing for <0 and exponentially for [0,1). 
%        .batchsize: 64 (defualt) Batchsize for neural network updates.
%        .buffersize: 1e6 (default) Experience buffer size
%        .noise: .1 (defualt) Perturbation noise radius for adv. and
%               set-based training methods.
%        .earlyStop: inf (default) Number of episodes for which the reward
%               did not increase before early stopping training.
%        .printFreq: 50 (default) Priniting frequency of training log.
%        .visRate: 50 (default) Visualisation frequency of learning
%               progress.
%        .actor.nn - Evaluation paramteres for the actor network
%           .poly_method: 'bounds' (default) Regression polynomial
%           .use_approx_error: true (default) Use approximation errors
%           .train - Training parameters for the actor network:
%               .use_gpu: true if available (default) Use CPU for training.
%               .optim: nnAdamOptimizer(1e-3,.9,.999,1e-8,1e-2) (default)
%                       Actor optimizer.
%               .backprop: true (default) Training boolean for actor.
%               .method: 'point'(default) Training method for actor:
%                   'point' Standard point-based training 
%                   'set' Set-based training [1]
%                   'random' Random adv. samples form perturbation ball
%                   'extreme' Adv. samples from edges of perturbation ball
%                   'naive' Adv. samples from naive algorithm [2]
%                   'grad' Adv. samples from grad algorthm [2]
%               .eta: 0.01 (default) Weighting factor for set-based 
%                   training of the actor.
%               .advOps - Parameters for adverserial training algs:
%                   .numSamples: 200 (default) Number of samples
%                   .alpha: 4 (default) Distribution paramater
%                   .beta: 4 (default) Distribution paramater
%               .zonotope_weight_update: 'outer_product' (default)
%                   Zonontope weight update for learnable params \theta
%        .critic.nn - Evaluation paramteres for the critic network:
%           .poly_method: 'bounds' (default) Regression polynomial
%           .use_approx_error: true (default) Use approximation errors
%           .train - Training parameters for the critic network:
%               .use_gpu: true if available (default) Use CPU for training.
%               .optim: nnAdamOptimizer(1e-3,.9,.999,1e-8) (default)
%                       Critic optimizer.
%               .backprop: true (default) Training boolean for critic.
%               .method: 'point'(default) Training method for critic:
%                   'point' Standard point-based training 
%                   'set' Set-based training [1]
%               .eta: 0.01 (default) Weighting factor for set-based 
%                   training of the critic.
%               .zonotope_weight_update: 'outer_product' (default)
%                   Zonontope weight update for learnable params \theta
%
% Outputs:
%   obj - genrated object
% 
% Refernces:
%   [1] Wendl, M. et al. Training Verifiably Robust Agents Using Set-Based 
%       Reinforcement Learning, 2024
%   [2] Pattanaik, A. et al. Robust Deep Reinforcement Learning with 
%       Adversarial Attacks, Int. Conf. on Autonomous Agents and Multiagent 
%       Systems (AAMAS) 2018
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: actor

% Authors:       Manuel Wendl
% Written:       24-October-2023
% Last update:   18-August-2024 (abstract class agentRL)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


properties
    
end

methods
    % contructor
    function obj = agentDDPG(actorNN,criticNN,varargin)
        % Parse function arguments
        narginchk(3,4)
        options = setDefaultValues(struct(),varargin);
        % init abstract class agentRL
        obj@agentRL(actorNN,criticNN,options);
    end

    function sobj = saveobj(obj)
        sobj = saveobj@agentRL(obj);
    end
end

methods (Access=protected)
    % trainNetworkStep - network training step
    function [obj,learnHistory] = trainNetworkStep(obj,randBatch, noiseBatchG, learnHistory, episode)
        % calculate critic target
        nextTargetAction = obj.targetActor.nn.evaluate_(randBatch{4},obj.options.rl.actor,obj.targetActor.idxLayer);
        targetQ = obj.targetCritic.nn.evaluate_([randBatch{4};nextTargetAction],obj.options.rl.critic,obj.targetCritic.idxLayer);
        targetBatch = randBatch{3} + (~randBatch{5}).*obj.options.rl.gamma.*targetQ;

        % learn critic
        [obj.critic,criticLoss] = obj.critic.train(randBatch,targetBatch,obj.options,noiseBatchG);

        learnHistory.criticLoss.center(episode) = learnHistory.criticLoss.center(episode) + criticLoss.center;
        learnHistory.criticLoss.vol(episode) = learnHistory.criticLoss.vol(episode) + criticLoss.vol;

        % learn actor
        [obj.actor, actorLoss] = obj.actor.train(obj.critic,randBatch,obj.options,noiseBatchG);

        learnHistory.actorLoss.center(episode) = learnHistory.actorLoss.center(episode) + actorLoss.center;
        learnHistory.actorLoss.vol(episode) = learnHistory.actorLoss.vol(episode) + actorLoss.vol;

        obj.targetCritic = obj.targetCritic.updateTargetCritic(obj.critic,obj.options);
        obj.targetActor = obj.targetActor.updateTargetActor(obj.actor,obj.options);
    end

    % (Potentially) gather weights of networks from GPU
    function obj = gatherNetworks(obj)
        obj.actor.nn.castWeights(single(1))
        obj.targetActor.nn.castWeights(single(1))
        obj.critic.nn.castWeights(single(1))
        obj.targetCritic.nn.castWeights(single(1))
    end

    % Delete gradients
    function obj = deleteAllGradients(obj)
        obj.actor.optim.deleteGrad(obj.actor.nn,obj.options);
        obj.targetActor.optim.deleteGrad(obj.targetActor.nn,obj.options);
        obj.critic.optim.deleteGrad(obj.critic.nn,obj.options);
        obj.targetCritic.optim.deleteGrad(obj.targetCritic.nn,obj.options);
    end
end
end

% ------------------------------ END OF CODE ------------------------------
