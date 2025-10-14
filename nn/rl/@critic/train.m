function [obj, loss] = train(obj,batch,targetBatch,options,noiseBatchG)
% train - train critc neural network
%   The actor is trained with randomly sampled batches from the experience
%   buffer. 
%
% Syntax:
%   [obj, loss] = train(obj,batch,targetBatch,options,noiseBatchG)
%
% Inputs:
%   batch - random batch from replay buffer 
%   targetBatch - targets of the batch
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
%               .eta: 0.01 (default) Weighting factor for set-based 
%                   training of the actor.
%               .advOps - Parameters for adverserial training algs:
%                   .numSamples: 200 (default) Number of samples
%                   .alpha: 4 (default) Distribution paramater
%                   .beta: 4 (default) Distribution paramater
%               .zonotope_weight_update: 'outer_product' (default)
%                   Zonontope weight update for learnable params \theta%   
%   noiseBatchG - pre-allocated noise batch
% 
% Outputs:
%   obj - updated critic
%   loss - mse loss of critic 
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
% Last update:   19-September-2024 (TL, renamed to train)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if strcmp(options.rl.critic.nn.train.method,'set')
    cBatch = cat(1,batch{1},reshape(batch{2}(:,1,:),size(batch{2}(:,1,:),1,3)));
    GBatch = cat(1,cat(2,noiseBatchG,zeros(size(batch{1},1),options.rl.critic.nn.train.num_init_gens-size(noiseBatchG,2),options.rl.batchsize,'like',noiseBatchG)),batch{2}(:,2:end,:));
else
    cBatch = cat(1,batch{1},batch{2});
    GBatch = [];
end

if ~strcmp(options.rl.critic.nn.train.method,'set')
    % Point based learning
    yPredBatch = obj.nn.evaluate_(cBatch,options.rl.critic,obj.idxLayer);
    [loss,gradOut,~] = aux_computeLoss(targetBatch,yPredBatch,[],options);
    obj.nn.backprop(gradOut,options.rl.critic,obj.idxLayer);
else
    % Set based Learinng
    [yPredBatchC,yPredBatchG] = obj.nn.evaluateZonotopeBatch_(cBatch,GBatch,options.rl.critic,obj.idxLayer);
    [loss,gradOutC,gradOutG] = aux_computeLoss(targetBatch,yPredBatchC,yPredBatchG,options);
    gradOutG = options.rl.critic.nn.train.eta * gradOutG;
    obj.nn.backpropZonotopeBatch_(gradOutC,gradOutG,options.rl.critic,obj.idxLayer,true);
end

obj.optim = obj.optim.step(obj.nn,options.rl.critic,obj.idxLayer);
end


% Auxiliary functions -----------------------------------------------------

function [loss,gradOutC,gradOutG] = aux_computeLoss(y,yPredC,yPredG,options)
% regression loss; half-squared error
lossFun = @(yPred,y) 1/(size(y,2)*2)*sum((y - yPred).^2,'all');
lossFunDer = @(yPred,y) - 1/size(y,2)*(y - yPred);
% compute loss of center
loss.center = lossFun(yPredC,y); 
% compute gradient of center
gradOutC = lossFunDer(yPredC,y);

if isempty(yPredG) % trained point-based
    loss.vol = 0; % no volume loss
    gradOutG = [];
else
    loss.vol = 1/size(yPredG,3)*sum(log(2*sum(abs(yPredG),2)),'all');

    gradOutG = 1/size(yPredG,3)*1/(sum(abs(yPredG),2)).*sign(yPredG);
    nanIdx = isnan(gradOutG)|isinf(gradOutG);
    gradOutG(nanIdx) = 0;
end

loss.vol =  1/max(options.rl.noise,[],'all')*loss.vol;
gradOutG = 1/max(options.rl.noise,[],'all')*gradOutG;
end

% ------------------------------ END OF CODE ------------------------------
