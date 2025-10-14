function [loss,grad_in] = getPolicyGradient(obj,batch,actionBatchC,actionBatchG,options,noiseBatchG)
% getPolicyGradient - get the policy gradient of the actor 
%
% Syntax:
%   [loss,grad_in] = getPolicyGradient(obj,batch,actionBatchC,actionBatchG,options,noiseBatchG)
%
% Inputs:
%   batch - random batch from replay buffer 
%   actionBatchC - centers of action batch
%   actionBatchG - generators of action batch
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
%   noiseBatchG - pre-allocated noise batch
%
% Outputs:
%   loss - actor policy gradient loss
%   grad_in - policy gradient and volume loss gradient of actor
% 
% Refernces:
%   [1] Wendl, M. et al. Training Verifiably Robust Agents Using Set-Based 
%       Reinforcement Learning, 2024
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: preprocessCriticInput
%
% See also: critic

% Authors:       Manuel Wendl
% Written:       03-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if all(strcmp('set',{options.rl.actor.nn.train.method,options.rl.critic.nn.train.method}))
    cBatch = cat(1,batch{1},actionBatchC);
    GBatch = cat(1,cat(2,noiseBatchG,zeros(size(batch{1},1),options.rl.critic.nn.train.num_init_gens-size(noiseBatchG,2),options.rl.batchsize,'like',noiseBatchG)),actionBatchG);
else
    cBatch = cat(1,batch{1},actionBatchC);
    GBatch = [];
end

% Get policy gtradient
if strcmp(options.rl.critic.nn.train.method,'point')
    % Point based learning
    yPredBatch = obj.nn.evaluate_(cBatch,options.rl.critic,obj.idxLayer);
    [loss, policyGradient, ~] = aux_computeLoss(yPredBatch,[]);
    grad_inC = obj.nn.backprop(policyGradient,options.rl.critic,obj.idxLayer,false);
    grad_inG = [];
else
    % Set based learning
    [yPredBatchC,yPredBatchG] = obj.nn.evaluateZonotopeBatch_(cBatch,GBatch,options.rl.critic,obj.idxLayer);
    [loss, policyGradient, gradOutG] = aux_computeLoss(yPredBatchC, yPredBatchG);
    [grad_inC,grad_inG] = obj.nn.backpropZonotopeBatch_(policyGradient,gradOutG,options.rl.critic,obj.idxLayer,false);
end

grad_in = aux_getActionInputGradients(grad_inC,grad_inG,batch,options);
end


% Auxiliary functions -----------------------------------------------------

function [loss, policyGradient, gradOutG] = aux_computeLoss(yPredC, yPredG)
% maximising the mean reward
lossFun = @(yPredC) - 1/size(yPredC,2) * sum(yPredC);
lossFunDer = @(yPredC) - 1/size(yPredC,2) * ones(size(yPredC));

% compute policy Gradient loss
loss.center = lossFun(yPredC);
% compute policy gradient 
policyGradient = lossFunDer(yPredC);

if isempty(yPredG)
    % Training point based.
    loss.vol = 0;
    gradOutG = [];
else
    % Compute set-based loss and gradient set.
    loss.vol = 1/size(yPredG,3)*sum(log(2*sum(abs(yPredG),2)),'all');
    gradOutG = 1/size(yPredG,3)*1/(sum(abs(yPredG),2)).*sign(yPredG);
    nanIdx = isnan(gradOutG)|isinf(gradOutG);
    gradOutG(nanIdx) = 0;
end

end

function grads_in_action = aux_getActionInputGradients(grad_inC,grad_inG,batch,options)
    grads_in_action.gradC = grad_inC(size(batch{1},1)+1:end,:);
    if ~isempty(grad_inG)
        grads_in_action.gradG = grad_inG(size(batch{1},1)+1:end,1:options.rl.critic.nn.train.num_init_gens,:);
    end
end

% ------------------------------ END OF CODE ------------------------------
