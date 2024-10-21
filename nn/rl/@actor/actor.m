classdef actor
% actor - actor of RLagent
%   The actor is given by a neural network of class neuralnetwork
%
% Syntax:
%   obj = actor(nn,options)
%
% Inputs:
%   nn - neuralNetwork of actor
%   options.rl - options for RL:
%        .noise: .1 (default) Perturbation noise radius for adv. and
%               set-based training methods.
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
%
% Outputs:
%   obj - generated actor
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
% Written:       23-October-2023
% Last update:   18-August-2024 (new allocation method)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    properties
        nn
        optim
        batchG
        idxLayer
    end

    methods
        % contructor
        function [obj, numGen, options] = actor(nn,options)
            % Parse function arguments
            narginchk(2,2)
            inputArgsCheck({ ...
                {nn,'att',{'neuralNetwork','neuralNetworkOld'}}; ...
                {options,'att','struct'}; ...
                })
            
            obj.nn = nn;

            % Obtain number of input dimension
            v0 = nn.neurons_in;
            % Validate input generators 
            if numel(options.rl.noise) == 1 % Case for l_inf ball
                options.rl.actor.nn.train.num_init_gens = v0;
            else % Case for custom noise
                options.rl.actor.nn.train.num_init_gens = size(options.rl.noise,2);
            end

            % layers
            obj.idxLayer = 1:length(obj.nn.layers);

            % Extract optimizer
            obj.optim = options.rl.actor.nn.train.optim;

            % To speed up computations and reduce gpu memory, we only use single
            % precision.
            inputDataClass = single(1);
            % Check if a gpu is used during training.
            useGpu = options.rl.actor.nn.train.use_gpu;
            if useGpu
                % Training data is also moved to gpu.
                inputDataClass = gpuArray(inputDataClass);
            end
            % (potentially) move weights of the network to gpu
            obj.nn.castWeights(inputDataClass);

            obj.optim.deleteGrad(obj.nn, options.rl.actor);

            % Allocate gpu memory: preallocate batch of generator matrices
            % This makes set-based training fast!
            if strcmp(options.rl.actor.nn.train.method,'set')
                % In each layer, store ids of active generators and identity matrices
                % for fast adding of approximation errors.
                numGen = obj.nn.prepareForZonoBatchEval(ones(v0,1),options.rl.actor,obj.idxLayer);
            else
                numGen = 0;
            end

        end

        function obj = updateTargetActor(obj,actor,options)
            % layers
            for i = obj.idxLayer
                if isa(obj.nn.layers{i},'nnLinearLayer')
                    % Update weigts and biases of target Actor
                    obj.nn.layers{i}.W = options.rl.tau*actor.nn.layers{i}.W + (1-options.rl.tau)*obj.nn.layers{i}.W;
                    obj.nn.layers{i}.b = options.rl.tau*actor.nn.layers{i}.b + (1-options.rl.tau)*obj.nn.layers{i}.b;
                end
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
