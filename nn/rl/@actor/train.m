function [obj, loss] = train(obj,critic,batch,options,noiseBatchG)
% train - actor neural network 
%   The actor is trained with randomly sampled batches from the experience
%   buffer. 
%
% Syntax:
%   [obj, loss] = train(obj,critic,batch,options,noiseBatchG)
%
% Inputs:
%   critic - critic with current parameters
%   batch - random batch from replay buffer
%   options.rl - options for RL:
%        .actor.nn - Evaluation paramteres for the actor network
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
%               .zonotope_weight_update: 'outer_product' (default)
%                   Zonontope weight update for learnable params \theta
%   noiseBatchG - Pre-allocated noise batch
%
% Outputs:
%   obj - updated actor
%   loss - actor loss 
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
% Written:       03-November-2023
% Last update:   19-September-2024 (TL, renamed to train)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
cBatch = batch{1};
GBatch = noiseBatchG;

if ~strcmp(options.rl.actor.nn.train.method,'set')
    % Point based learning
    actionBatchC = obj.nn.evaluate_(cBatch,options.rl.actor,obj.idxLayer);
    [loss,policyGradient] = critic.getPolicyGradient(batch,actionBatchC,[],options);
    obj.nn.backprop(policyGradient.gradC,options.rl.actor,obj.idxLayer);
    loss.vol = 0;
else
    % Set based Learinng
    [actionBatchC,actionBatchG] = obj.nn.evaluateZonotopeBatch_(cBatch,GBatch,options.rl.actor,obj.idxLayer);
   
    [loss,policyGradient] = critic.getPolicyGradient(batch,actionBatchC,actionBatchG,options,noiseBatchG);

    if strcmp(options.rl.critic.nn.train.method,'point') 
        [loss.vol,gradOutG] = aux_computeVolumeLoss(actionBatchG);
    elseif strcmp(options.rl.critic.nn.train.method,'set')
        [loss.vol,gradOutG] = aux_computeVolumeLoss(actionBatchG);
        gradOutG = options.rl.actor.nn.train.omega * gradOutG + (1-options.rl.actor.nn.train.omega) * policyGradient.gradG;
    else
        throw(CORAerror("CORA:notDefined",'Other trainig methods than point or set are not implemented for the critic.'))
    end

    % Scale volume gradients
    gradOutG = options.rl.actor.nn.train.eta/max(options.rl.noise,[],'all') * gradOutG;
    loss.vol =  1/max(options.rl.noise,[],'all')*loss.vol;

    obj.nn.backpropZonotopeBatch_(policyGradient.gradC,gradOutG,options.rl.actor,obj.idxLayer);

end

obj.optim = obj.optim.step(obj.nn,options.rl.actor,obj.idxLayer);

end


% Auxiliary functions -----------------------------------------------------

function [loss,gradOutG] = aux_computeVolumeLoss(yPredG)
loss = 1/size(yPredG,3)*sum(log(2*sum(abs(yPredG),2)),'all');
gradOutG = 1/size(yPredG,3)*1/(sum(abs(yPredG),2)).*sign(yPredG);
nanIdx = isnan(gradOutG)|isinf(gradOutG);
gradOutG(nanIdx) = 0;
end

% ------------------------------ END OF CODE ------------------------------
