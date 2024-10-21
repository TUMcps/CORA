classdef critic
% critic - critic of RLagent
%   The critic is given by a neural network of class neuralnetwork
%
% Syntax:
%   obj = critic(nn,numGensActor,options)
%
% Inputs:
%   nn - neuralNetwork of critic
%   numGensActor - number of generators of the actor output
%   options.rl - options for RL:
%        .noise: .1 (default) Perturbation noise radius for adv. and
%               set-based training methods.
%        .critic.nn - Evaluation paramteres for the critic network
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
%
% Outputs:
%   obj - generated critic
%
% Refernces:
%   [1] Wendl, M. et al. Training Verifiably Robust Agents Using Set-Based 
%       Reinforcement Learning, 2024
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: DDPGagent

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
        function [obj, options] = critic(nn,numGensActor,options) 
            % Parse function arguments
            narginchk(3,3)            
            inputArgsCheck({ ...
                {nn,'att',{'neuralNetwork','neuralNetworkOld'}}; ...
                {numGensActor, 'att', 'numeric'}; ...
                {options,'att','struct'}; ...
                })

            obj.nn = nn;

            % Obtain number of input dimension
            v0 = nn.neurons_in;
            % Set number of input generators
            options.rl.critic.nn.train.num_init_gens = numGensActor;

            % layers
            obj.idxLayer = 1:length(obj.nn.layers);

            % Extract optimizer
            obj.optim = options.rl.critic.nn.train.optim;

            % To speed up computations and reduce gpu memory, we only use single
            % precision.
            inputDataClass = single(1);
            % Check if a gpu is used during training.
            useGpu = options.rl.critic.nn.train.use_gpu;
            if useGpu
                % Training data is also moved to gpu.
                inputDataClass = gpuArray(inputDataClass);
            end
            % (potentially) move weights of the network to gpu
            obj.nn.castWeights(inputDataClass);

            obj.optim.deleteGrad(obj.nn,options.rl.critic);

            % Allocate gpu memory: preallocate batch of generator matrices
            % This makes set-based training fast!
            if strcmp(options.rl.critic.nn.train.method,'set')
                % In each layer, store ids of active generators and identity matrices
                % for fast adding of approximation errors.
                numGen = obj.nn.prepareForZonoBatchEval(ones(v0,1),options.rl.critic,obj.idxLayer);
                % Allocate generators for initial perturbance set.
                idMat = eye(v0,'like',inputDataClass);
                obj.batchG = cast(repmat(idMat,1,1,options.rl.batchsize),'like',inputDataClass);
            end
        end

        function obj = updateTargetCritic(obj,critic,options)
            % layers
            for i = obj.idxLayer
                if isa(obj.nn.layers{i},'nnLinearLayer')
                    % Update weigts and biases of target Critic
                    obj.nn.layers{i}.W = options.rl.tau*critic.nn.layers{i}.W + (1-options.rl.tau)*obj.nn.layers{i}.W;
                    obj.nn.layers{i}.b = options.rl.tau*critic.nn.layers{i}.b + (1-options.rl.tau)*obj.nn.layers{i}.b;
                end
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
