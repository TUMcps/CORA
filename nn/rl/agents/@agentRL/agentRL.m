classdef (Abstract) agentRL
% agentRL - abstract agentRL
%
% Syntax:
%   obj = agentRL(actorNN,criticNN,options)
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
%   obj - generated object
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
% Written:       18-August-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    actor
    critic
    targetActor
    targetCritic
    buffer
    options
    learnHistory
end

methods
    % contructor
    function obj = agentRL(actorNN,criticNN,options)
        % Validate function arguments
        inputArgsCheck({ ...
            {actorNN,'att','neuralNetwork'}; ...
            {criticNN,'att','neuralNetwork'}; ...
            {options,'att','struct'}; ...
            });
        % Set default training parameters
        options = nnHelper.validateRLoptions(options);

        % instantiate actor and critic
        [obj.actor, numGensActor, options] = actor(actorNN.copyNeuralNetwork(),options);
        [obj.critic, options] = critic(criticNN.copyNeuralNetwork(),numGensActor,options);

        % Instantiate target networks
        [obj.targetActor, numGensActor, options] = actor(actorNN.copyNeuralNetwork(),options);
        [obj.targetCritic, options] = critic(criticNN.copyNeuralNetwork(),numGensActor,options);

        % instantiate buffer
        obj.buffer = buffer(options.rl.buffersize);

        % store options
        obj.options = options;
    end

    function sobj = saveobj(obj)
        sobj = obj;
        sobj.buffer.array = [];
        s = rng;
        sobj.options.rl.trainSeed = s.Seed;
    end

end

methods  (Access=protected, Abstract)
    % Training
    [obj,learnHistory] = trainNetworkStep(obj,randBatch, noiseBatchG, learnHistory, episode)
    % Remove networks from the GPU
    obj = gatherNetworks(obj)
    % Clear gradients
    obj = deleteAllGradients(obj)
end


end

% Auxiliary functions -----------------------------------------------------

% ------------------------------ END OF CODE ------------------------------
