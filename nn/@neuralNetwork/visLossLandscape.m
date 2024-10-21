function [alpha,beta,lossCenter,lossVol,lossTotal] = ...
    visLossLandscape(obj,x,t,options,varargin)
% visLossLandscape - Visualize the loss landscape of a neural network.
%
% Syntax:
%    [alpha,beta,lossCenter,lossVol,lossTotal] =...
%       nn.visLossLandscape(x,t,options,numSamples);
%    figure;
%    contour(alpha,beta,lossCenter,10); % contour plot of loss landscape
%
% Inputs:
%    obj - neural network
%    x - test inputs; i-th input x(:,i)
%    t - test targets; i-th target t(:,i)
%    options - training parameters (see neuralNetwork.train),
%           needed to pass information how loss is computed (training
%           perturbation radius, scaling of volume in loss, etc.)
%    numSamples - number of samples along gradient axis (default=20)
%
% Outputs:
%    alpha - values along the first axis; represents a gradient direction 
%    beta - values along the second axis; represents a gradient direction 
%    lossCenter - matrix with loss values, center error
%    lossVol - matrix with loss values, volume heuristic error
%    lossTotal - matrix with loss values, total error
%    
% References:
%    [1] Hao Li, Zheng Xu, Gavin Taylor, Christoph Studer, et al. 
%        "Visualizing the loss landscape of neural nets". In NeurIPS, 
%        pages 6391–6401, 2018. 4, 5
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Lukas Koller
% Written:       21-June-2023
% Last update:   02-August-2023 (adapted variable names)
% Last revision: ---    

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(4,5)
[numSamples] = setDefaultValues({20},varargin);

% validate input
inputArgsCheck({ ...
    {obj,'att','neuralNetwork'}; ...
    {x,'att','numeric','array'}; ...
    {t,'att','numeric','array'}; ... 
    {options,'att','struct'}; ... 
    {numSamples,'att','numeric'};
})

% use the neuralNetwork.train function to compute the loss
% setup training params
options.nn.train.optim = nnSGDOptimizer(0);
options.nn.train.max_epoch = 1;

% select data samples
% idx = randsample(size(x,2),trainParamsVis.mini_batch_size);
% x = x(:,idx);
% t = t(:,idx);

rng('default')
% sample random weight directions \delta and \eta
delta = aux_generateRandomNormalizedWeights(obj);
eta = aux_generateRandomNormalizedWeights(obj);
% combine sampled weight direction with volume gradient
% zeta = aux_generateNormalizedLossGradient(obj,x,t,trainParamsVis);
% delta = aux_convexCombOfNetworkParams(delta,zeta,0.1);
% eta = aux_convexCombOfNetworkParams(eta,zeta,0.1);

% scales for plot
alpha = linspace(-1.5,1.5,numSamples);
beta = linspace(-1.5,1.5,numSamples);

% store loss
lossTotal = zeros(length(alpha),length(beta));
lossCenter = zeros(length(alpha),length(beta));
lossVol = zeros(length(alpha),length(beta));
% compute losses
for i=1:length(alpha)
    for j=1:length(beta)
        % create copy of the neural network to modify the weights
        nnTest = copyNeuralNetwork(obj);
        for k=1:length(nnTest.layers)
            % extract original weights and biases
            layerk = obj.layers{k};
            if isa(layerk,'nnLinearLayer') | isa(layerk,'nnConv2dLayer')
                Wk = layerk.W;
                bk = layerk.b;
                % compute weight and bias updates
                deltak = delta{k};
                etak = eta{k};
                dWk = alpha(i)*deltak(:,1:end-1) + beta(j)*etak(:,1:end-1);
                dbk = alpha(i)*deltak(:,end) + beta(j)*etak(:,end);
                % update weight matrix and bias
                nnTest.layers{k}.W = Wk + dWk;
                nnTest.layers{k}.b = bk + dbk;
            end
        end
        % compute loss for input data (just validation data)
        lossij = nnTest.train(x,t,[],[],options,false);
        lossCenter(i,j) = lossij.center; % lossij.valCenter;
        lossVol(i,j) = lossij.vol; % lossij.valVol;
        lossTotal(i,j) = lossij.total; % lossij.valTotal;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function delta = aux_generateRandomNormalizedWeights(nn)
% copy neural network
nnCopy = copyNeuralNetwork(nn);
% initialize random weights
% nnCopy.initWeights('glorot');

numLayers = length(nnCopy.layers);
% sample random weight directions \delta and \eta
delta = cell(numLayers,1);
for i=1:numLayers
    layeri = nn.layers{i};
    if isa(layeri,'nnLinearLayer') | isa(layeri,'nnConv2dLayer')
        % extract weight matrix and bias
        Wi = layeri.W;
        bi = layeri.b;
        [m,n] = size(Wi);
        % sample random directions with appropriate size + normalize
        deltaWi = mvnrnd(zeros(m,1),eye(m),n)';
        deltaWi = norm(Wi,'fro')/norm(deltaWi,'fro')*deltaWi; % normalize
        deltabi = mvnrnd(zeros(m,1),eye(m),1)';
        deltabi = norm(bi,'fro')/norm(deltabi,'fro')*deltabi; % normalize
        delta{i} = [deltaWi deltabi];
        % delta{i} = [Wi bi];
    else
        delta{i} = [];
    end
end
end

function delta = aux_generateNormalizedLossGradient(nn,x,t,trainParams)
[v0,~] = size(x);

numInitGens = trainParams.num_init_gens;
numApproxErr = trainParams.num_approx_err;
% copy neural network
nnCopy = copyNeuralNetwork(nn);
% initialize network for batch propagation
nnCopy.prepareForZonoBatchEval(x,min(v0,numInitGens),numApproxErr);

% init output
numLayers = length(nnCopy.layers);
delta = cell(numLayers,1);

% propagate inputs through the network
nnCopy.train(x,t,[],[],trainParams,false);

% extract gradient from network
for i=1:numLayers
    layeri = nnCopy.layers{i};
    if isa(layeri,'nnLinearLayer') | isa(layeri,'nnConv2dLayer')
        % extract weight matrix and bias
        Wi = layeri.W;
        bi = layeri.b;
        % extract gradient for weight matrix and bias
        deltaWi = layeri.backprop.grad.W;
        deltaWi = norm(Wi,'fro')/norm(deltaWi,'fro')*deltaWi; % normalize
        deltabi = layeri.backprop.grad.b;
        deltabi = norm(bi,'fro')/norm(deltabi,'fro')*deltabi; % normalize
        delta{i} = [deltaWi deltabi];
    else
        delta{i} = [];
    end
end
end

function delta = aux_convexCombOfNetworkParams(delta1,delta2,lambda)
numLayers = length(delta1);
delta = cell(numLayers,1);
for i=1:numLayers
    delta{i} = lambda*delta1{i} + (1-lambda)*delta2{i};
end
end

% ------------------------------ END OF CODE ------------------------------
