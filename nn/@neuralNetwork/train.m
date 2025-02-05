function [loss,trainTime] = train(nn, trainX, trainT, valX, valT, varargin)
% train - train a neural network
%
% Syntax:
%    % initialize the weights of the neural network
%    nn.initWeights('glorot');
%    % specify training parameters
%    options.nn.train = struct(...
%       'use_gpu',true, 'optim',nnAdamOptimizer(1e-3), 'method','point');
%    % train the neural networks
%    [loss,trainTime] = nn.train(trainX, trainT, valX, valT, options);
%
% Inputs:
%    nn - neural network
%    trainX - training input; i-th input trainX(:,i)
%    trainT - training targets; i-th target trainT(:,i)
%    valX - validation input, normalized to [0,1]
%    valT - validation targets
%    options.nn.train - possible training parameters
%        .optim: instance of nnOptimizer, e.g. nnSGDOptimizer,
%           nnAdamOptimizer
%        .max_epoch: maximum number of training epochs
%        .mini_batch_size: training batch size
%        .loss: type of loss function, {'mse','softmax+log','custom'};
%           'mse' for regression loss, 'softmax+log' for cross-entropy, 
%               or 'custom': in this case provide .loss_fun and .loss_grad
%        .noise: training perturbatuion noise (default: 0.0)
%        .input_space_inf: lower bound of the input space (default: 0)
%        .input_space_sup: upper bound of the input space (default: 1)
%        .warm_up: number of training epochs at the beginning without any 
%           perturbation noise, (default: 0)
%        .ramp_up: linearly ramp up perturbation noise from warm-up to
%           ramp-up, (default: 0)
%        .method: training methods, {'point','set','madry','gowal','trades','sabr'}; 
%           extra options for the individual
%           'point': regular point-based training, see [2]
%           'set': set-based training
%               .volume_heuristic: {'interval','f-radius'}
%               .tau: scaling-factor of volume heuristic in loss
%               .zonotope_weight_update: {'center','sample','extreme','sum'}
%               .num_approx_err: number of approx. errors per activation 
%                   layer during training, inf for all
%               .num_init_gens: number of generators of the input 
%                   zonotopes; ignored if initial_generators='l_inf'
%               .init_gens: {'l_inf','random','sensitivity'}
%               .poly_method: {'bounds','singh','center'}
%           'madry': augment training inputs with adversarial attacks, see [3]
%               .pgd_iterations: number of iterations
%               .pgd_stepsize: perturbation stepsize (default: 0.01)
%           'gowal': interval-bound propagation, see [4]
%               .kappa: weighting factor (final)
%           'trades': tradeoff (betw. accuracy & robustness) inspired
%               training, see [6]
%               .lambda: weighting factor
%               .pgd_iterations: number of iterations
%               .pgd_stepsize: perturbation stepsize (default: 0.01)
%           'sabr': combination of adversarial attacks and interval-bound
%               propagation, see [5]
%               .lambda: weighting factor
%               .pgd_iterations: number of iterations
%               .pgd_stepsize: perturbation stepsize
%               .pgd_stepsize_decay: decay factor of stepsize in PGD
%               .pgd_stepsize_decay_iter: decay iterations of stepsize in PGD
%        .shuffle_data: shuffle training data, {'never','every_epoch'},
%           (default: 'never')
%        .lr_decay: decay of learning rate (default: 1.0)
%        .lr_decay_epoch: epoch at which learning rate is decayed
%           (default: [])
%        .early_stop: number, abort training if validation loss is 
%           non-decreasing after certain amount of steps (default: inf)
%        .val_freq: validation frequency (default: 50)
%        .print_freq: print losses every printFreq-th iteration 
%           (default: 50)
%    verbose - print intermediate training losses (default: true)
%
% Outputs:
%    nn - trained neural network
%    loss - loss during training, struct
%        .center: (training) loss
%        .vol: (training) volume heuristic loss
%        .total: total training training loss
%        .val: validation loss
%    trainTime - training time
%    
% References:
%    [1] C. M. Bishop. Pattern Recognition and Machine Learning. 
%        Information Science and Statistics. Springer New York, NY, 2006.
%    [2] Koller, L. "Co-Design for Training and Verifying Neural Networks",
%           Master's Thesis
%    [3] Madry, A. et al. Towards deep learning models resistant to 
%           adversarial attacks. ICLR. 2018
%    [4] Gowal, S. et al. Scalable verified training for provably
%           robust image classification. ICCV. 2019
%    [5] Mueller, M. et al. Certified Training: Small Boxes are All You 
%           Need. ICLR. 2022
%    [6] Zhang, H. et al. Theoretically Principled Trade-off between 
%           Robustness and Accuracy. ICML. 2019
% 
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork, neuralNetwork/initWeights

% Authors:       Lukas Koller
% Written:       22-May-2023
% Last update:   05-June-2023 (modified signature)
%                05-July-2023 (clean up + add early-stopping)
%                25-July-2023 (zero gradients before training)
%                31-July-2023 (parallelization, point-wise batch-training, cross-entropy loss)
%                02-August-2023 (added batch-eval & -backprop for zonotope)
%                19-August-2023 (memory optimizations for GPU training)
%                16-January-2024 (sensitivty-based generators)
%                22-January-2024 (added IBP-based & TRADES training method)
%                07-February-2024 (added SABR training method + simplified f-radius implementation)
%                22-February-2024 (trainParams, merged options.nn)
% Last revision: ---    

% ------------------------------ BEGIN CODE -------------------------------

% Parse function arguments
narginchk(5,7)
[options,verbose] = setDefaultValues({struct,true},varargin);
% Validate function arguments
inputArgsCheck({ ...
    {nn,'att','neuralNetwork'}; ...
    {trainX,'att','numeric','array'}; ...
    {trainT,'att','numeric','array'}; ... 
    {valX,'att','numeric','array'}; ... 
    {valT,'att','numeric','array'}; ... 
    {options,'att','struct'}; ... 
    {verbose,'att','logical'};
});
% Set default training parameters
options = nnHelper.validateNNoptions(options,true);

% Extract training parameters
maxEpoch = options.nn.train.max_epoch;
miniBatchSize = options.nn.train.mini_batch_size;

maxNoise = options.nn.train.noise;
% Compute numbre of ramp up epochs
numRampUpEpochs = options.nn.train.ramp_up - options.nn.train.warm_up;

% Obtain number of input dimension and number of data points
[v0,N] = size(trainX);
% Validate input generators
if strcmp(options.nn.train.init_gens,'l_inf')
    options.nn.train.num_init_gens = v0;
else
    options.nn.train.num_init_gens = min(v0,options.nn.train.num_init_gens);
end

% enable backpropagation
options.nn.train.backprop = true;

% Compute indices of iterations where to validate and/or print verbose output
numIter = maxEpoch*max(1,floor(N/miniBatchSize)); % total number of iterations
iterPerEpoch = ceil(N/miniBatchSize); % number of iteration per epoch
iter = 1; % counter for current iteration
printFreq = options.nn.train.print_freq;
% Compute iterations at which the current loss is printed
printIter = unique([1,printFreq:printFreq:numIter,numIter]); 
% Compute iterations where a validate is done
valFreq = options.nn.train.val_freq;
if valFreq < 0 | isempty(valX)
    valIter = [];
else
    valIter = unique([1,valFreq:valFreq:numIter,numIter]); 
end

% Extract optimizer
optim = options.nn.train.optim;
optim.lrDecay = options.nn.train.lr_decay;
optim.lrDecayIter = iterPerEpoch*options.nn.train.lr_decay_epoch;


ivalC = options.nn.interval_center;

% Print training parameters
if verbose
    aux_printParameters(options,optim,maxEpoch,miniBatchSize,maxNoise)
end

% Initialize loss struct to store loss values
loss = struct('center',zeros(numIter,1),...
    'vol',zeros(numIter,1),...
    'total',zeros(numIter,1),...
    'valCenter',zeros(length(valIter),1),...
    'valVol',zeros(length(valIter),1),...
    'valTotal',zeros(length(valIter),1)...
);

% Count number of non-decreasing validation steps for early stopping
minValLoss = Inf;
numNonDecrVal = 0;

% Delete gradients
optim.deleteGrad(nn,options);

% Get indices of layers for propagation. The final softmax layer is
% considered as part of the loss function.
idxLayer = 1:length(nn.layers);
if strcmp(options.nn.train.loss,'softmax+log')
    idxLayer = idxLayer( ...
        arrayfun(@(i) ~isa(nn.layers{i},'nnSoftmaxLayer'),idxLayer));
end

% To speed up computations and reduce gpu memory, we only use single 
% precision.
inputDataClass = single(1);
% Check if a gpu is used during training.
useGpu = options.nn.train.use_gpu;
if useGpu
    % Training data is also moved to gpu.
    inputDataClass = gpuArray(inputDataClass);
end
% (potentially) move weights of the network to gpu
nn.castWeights(inputDataClass);

% Allocate gpu memory: preallocate batch of generator matrices
% This makes set-based training fast!
if strcmp(options.nn.train.method,'set')
    % In each layer, store ids of active generators and identity matrices 
    % for fast adding of approximation errors.
    numGen = nn.prepareForZonoBatchEval(trainX,options,idxLayer);
    % Allocate generators for initial perturbance set.
    idMat = eye(v0,'like',inputDataClass);
    batchG = cast(repmat(idMat,1,1,miniBatchSize),'like',inputDataClass);
end

% set up verbose output table
table = aux_setUpTrainingTable(options);
if verbose
    table.printHeader();
end

% Start training time
tic

% Main training loop
for epoch=1:maxEpoch % epoch
    if strcmp('every_epoch',options.nn.train.shuffle_data)
        % permute inputs
        pIdx = randperm(N);
        trainX = trainX(:,pIdx);
        trainT = trainT(:,pIdx);
    end

    for i=1:miniBatchSize:max(1,N - mod(N,miniBatchSize)) % training iterations
        % Compute a scaling factor to scale training noise based on warm-up
        % and ramp-up epochs.
        iterScale = (iter - options.nn.train.warm_up*iterPerEpoch)...
            /(numRampUpEpochs*iterPerEpoch);
        iterScale = min(max(0,iterScale),1);
        % Scale training noise
        options.nn.train.noise = maxNoise*iterScale;

        % Compute indices for elements in current batch
        miniBatchIdx = i:i+min(miniBatchSize,N)-1;
        % Extract data of current batch
        xBatch = cast(trainX(:,miniBatchIdx),'like',inputDataClass);
        tBatch = cast(trainT(:,miniBatchIdx),'like',inputDataClass);
        % enable backpropagation
        options.nn.train.backprop = true;
        options.nn.interval_center = ivalC;
        if strcmp(options.nn.train.method,'point')
            % Forward propagation
            yBatch = nn.evaluate_(xBatch,options,idxLayer);
            % Compute loss
            [lossBatch,grad,~] = aux_computeLoss(tBatch,yBatch,[],options);
            % Backpropagation
            nn.backprop(grad,options,idxLayer);
        elseif strcmp(options.nn.train.method,'set')
            % construct input generators
            [xBatch,xBatchG] = aux_construtInputGenerators(nn,xBatch, ...
                batchG(:,:,1:length(miniBatchIdx)),idMat,tBatch,options);
            % Set-based forward propagation
            [ycBatch,yGBatch] = nn.evaluateZonotopeBatch_(xBatch,xBatchG, ...
                options,idxLayer);
            % Compute set-based loss
            [lossBatch,gc,gG] = aux_computeLoss(tBatch,ycBatch,yGBatch, ...
                options);
            % Set-based backpropagation
            nn.backpropZonotopeBatch_(gc,gG,options,idxLayer);
        else
            % Train using other methods ('madry','gowal','trades','sabr'), 
            % see [3-6]
            lossBatch = aux_trainOtherMethods(nn,xBatch,tBatch,idxLayer,...
                options,iterScale);
        end
        % Store loss of current batch
        loss.center(iter) = lossBatch.center;
        loss.vol(iter) = lossBatch.vol;
        loss.total(iter) = lossBatch.total;
        % To save memory, we clear all variables that are only used
        % during the batch computation
        batchVars = {'xBatch','tBatch','yBatch','grad','xBatchG','ycBatch',...
            'yGBatch','gc','gG'};
        clear(batchVars{:});

        % Step the optimizer to update the weights
        optim = optim.step(nn,options);
        
        if ismember(iter,valIter) 
            % Run validation
            valN = size(valX,2);
            % compute number of validation batches
            numValBatches = ceil(valN/miniBatchSize);
            % initialize array for validation losses
            lossValCenter = zeros(numValBatches,1);
            lossValVol = zeros(numValBatches,1);
            lossValTotal = zeros(numValBatches,1);
            % disable backpropagation
            options.nn.train.backprop = false;
            options.nn.interval_center = false;
            for j=1:miniBatchSize:max(1,valN - mod(valN,miniBatchSize))
                % compute indices for elements in current batch
                valBatchIdx = j:min(j+miniBatchSize-1,valN);
                % extract data of current batch
                valXBatch = cast(valX(:,valBatchIdx),'like',inputDataClass);
                valTBatch = cast(valT(:,valBatchIdx),'like',inputDataClass);
                % Point-based forward propagation
                valYBatch = nn.evaluate_(valXBatch,options,idxLayer);
                % Compute loss
                [lossValBatch,~,~] = aux_computeLoss(valTBatch,valYBatch,[],options);
                % Store loss
                lossValCenter(j) = lossValBatch.center*length(valBatchIdx);
                lossValVol(j) = lossValBatch.vol*length(valBatchIdx);
                lossValTotal(j) = lossValBatch.total*length(valBatchIdx);

                % To save memory, we clear all variables that are only used
                % during the batch computation
                batchVars = {'valXBatch','valTBatch','valYBatch','lossValBatch'};
                clear(batchVars{:});
            end
            % Store validation loss
            valIterIdx = find(valIter == iter);
            loss.valCenter(valIterIdx) = sum(lossValCenter)/valN;
            loss.valVol(valIterIdx) = sum(lossValVol)/valN;
            loss.valTotal(valIterIdx) = sum(lossValTotal)/valN;
            % Early stopping: check if validation loss is decreasing 
            if minValLoss > loss.valTotal(valIterIdx)
                numNonDecrVal = 0;
                minValLoss = loss.valTotal(valIterIdx);
            else
                numNonDecrVal = numNonDecrVal + 1;
            end
        end

        % Print loss of current training iteration
        if verbose && ismember(iter,printIter) 
            trainTime = toc;
            aux_printTrainingIteration(table,iter,epoch,trainTime,valIter,loss)
        end

        % Check for early-stopping
        if numNonDecrVal > options.nn.train.early_stop
            if verbose
                fprintf(['Stopped early! The loss for the last %i ' ...
                    'validation steps did not decrese.\n'], ...
                        options.nn.train.early_stop);
            end
            break;
        end

        % Increment counter for current training iteration
        iter = iter + 1;
    end
    % Stop training time
    trainTime = toc;
    % Check if we stopped early
    if numNonDecrVal > options.nn.train.early_stop
        break;
    end
end

% brint bottom
if verbose
    table.printFooter();
end

% Clear variables that are no longer used.
vars = {'batchG','idMat'};
clear(vars{:});

% (potentially) gather weights of the network from gpu
nn.castWeights(single(1));

% Delete gradients
optim.deleteGrad(nn,options);

end


% Auxiliary functions -----------------------------------------------------

function [loss,gc,gG] = aux_computeLoss(t,yc,yG,options)
% Compute the loss for a given zonotope Y=(yc,yG). If the generator matrix
% is empty, we just compute the point based loss.

% Obtain mini-batch size
[vK,batchSize] = size(t);

% Retrieve loss function
if strcmp(options.nn.train.loss,'mse') 
    % Regression loss (half-squared error)
    lossFun = @(t,y) 1/batchSize*0.5*sum((y - t).^2,'all');
    gradFun = @(t,y) 1/batchSize*(y + (-t));
elseif strcmp(options.nn.train.loss,'softmax+log')
    % Improve numerical stability of softmax
    % yc = yc - max(yc,[],1);
    % Use numerically stable log-softmax
    logSoftmax = @(x) x - log(sum(exp(x),1));
    % Classification loss (cross-entropy)
    lossFun = @(t,y) -1/batchSize*sum(t.*logSoftmax(y),'all');
    gradFun = @(t,y) 1/batchSize*(softmax(y) + (-t));
elseif strcmp(options.nn.train.loss,'custom')
    % Take given loss function and derivative
    lossFun = options.nn.train.loss_fun;
    gradFun = options.nn.train.loss_grad;
else
    throw(CORAerror('CORA:wrongFieldValue', ...
        sprintf("Unsported loss function: %s; Only supported values " + ...
            "for 'point_loss_fun_type' are 'mse','softmax+log', and " + ...
                "'custom'!",options.nn.train.loss)));
end

if options.nn.interval_center
    ycl = reshape(yc(:,1,:),[vK batchSize]);
    ycu = reshape(yc(:,2,:),[vK batchSize]);

    % % Gowal-IBP Loss
    % z = t.*ycl + (1-t).*ycu;
    % % Compute loss of center
    % loss.center = cast(lossFun(t,z),'like',single(1));
    % % Compute gradient of center
    % gc = gradFun(t,z);
    % gc = permute(cat(3,t.*gc,(1-t).*gc),[1 3 2]);

    % Center Loss
    z = 1/2*(ycu + ycl);
    % Compute loss of center
    loss.center = cast(lossFun(t,z),'like',single(1));    
    % Compute gradient of center
    gc = gradFun(t,z);
    gc = 1/2*permute(cat(3,gc,gc),[1 3 2]);

    idMat = eye(size(z,1),'like',z);
    yG = cat(2,yG,1/2*permute(ycu - ycl,[1 3 2]).*idMat);
else
    % Use center.
    z = yc;
    % Compute loss of center.
    loss.center = cast(lossFun(t,z),'like',single(1));
    % Compute gradient of center.
    gc = gradFun(t,z);
end


if isempty(yG) 
    % Training point-based: There is no volume loss
    loss.vol = 0;
    loss.total = loss.center;
    % There is not gradient for the generator matrix
    gG = zeros(size(yG),'like',yG);
else 
    % Compute volume heuristic
    if strcmp(options.nn.train.volume_heuristic,'interval')
        % The interval norm is the sum of all absolute values of the
        % generator matrix.
        yPredVolHeu = 1/batchSize*sum(abs(yG),'all');
        gG = 1/batchSize*sign(yG);
    elseif strcmp(options.nn.train.volume_heuristic,'f-radius')
        % The F-radius is the square-root of the sum of all squared entries
        % of the generator matrix.
        frad = sqrt(sum(yG.^2,[1 2]));
        zIdx = (frad(:) < eps('like',frad));
        yPredVolHeu = 1/batchSize*sum(frad);
        gG = 1/batchSize*(yG./frad);
        gG(:,:,zIdx) = 0;
    else
        throw(CORAerror('CORA:wrongFieldValue', ...
            sprintf("Unsported volume heuristic: %s; Only supported values " + ...
                "for 'volume_heuristic' are 'interval' and 'f-radius'!!",...
                    options.nn.train.volume_heuristic)));
    end
    loss.vol = cast(yPredVolHeu,'like',single(1));
    % Extract training parameters
    noise = options.nn.train.noise;
    tau = options.nn.train.tau;
    % Normalize volume loss
    if noise ~= 0
        loss.vol = 1/(vK*noise)*loss.vol;
        gG = 1/(vK*noise)*gG;
    else
        loss.vol = 0;
        gG = zeros(size(gG),'like',gG);
    end
    % Compute weights to combine center and volume loss
    cCenter = 1 - tau;
    cVol = tau;
    % Compute total loss
    loss.total = cCenter*loss.center + cVol*loss.vol;
    % Scale the gradients according to their weights
    gc = cCenter*gc;
    gG = cVol*gG; 

    if strcmp(options.nn.train.zonotope_weight_update,'outer_product')
        % For 'outer_product' or 'sum', we compute the set of 
        % classification losses. Thus, we additionally add the generator
        % matrix to the gradient.
        gG = gG + 1/batchSize*cCenter*yG;
    elseif strcmp(options.nn.train.zonotope_weight_update,'sum')
        % Ensure uniform distribution over the output set.
        % yGnorm = yG./max(vecnorm(yG),eps('like',yG));
        % Ebbt = vK/(vK+2)*pagemtimes(yGnorm,'transpose',yGnorm,'none');
        % Ebbt = vK/(vK+2)*sign(pagemtimes(yG,'transpose',yG,'none'));
        % Ebbt = Ebbt./max(Ebbt,[],3);
        % Ebbt(isnan(Ebbt)) = 1;

        % gG = gG + pagemtimes(1/batchSize*cCenter*yG,Ebbt);
        % gG = pagemtimes(gG + 1/batchSize*cCenter*yG,Ebbt);
        % gG = gG + 1/batchSize*cCenter*yG;
    end

    if options.nn.interval_center
        q = size(yG,2) - vK;
        diagIdx = reshape(sub2ind(size(yG),repmat(1:vK,1,batchSize),...
            q + repmat(1:vK,1,batchSize),repelem(1:batchSize,1,vK)),...
            [vK batchSize]);
        gc = gc + 1/2*permute(cat(3,-gG(diagIdx),gG(diagIdx)),[1 3 2]);
        gG = gG(:,1:q,:); % 1/3.*yG(:,1:q,:); % 
    end
end

end

function [xBatch,xBatchG] = aux_construtInputGenerators(nn,xBatch,batchG,idMat,tBatch,options)
% Construct the generator matrix for the input zonotopes.
[v0,batchSize] = size(xBatch);
% Extract training parameters
epsilon = options.nn.train.noise;
% Obtain the bound of the input space.
inpInf = options.nn.train.input_space_inf;
inpSup = options.nn.train.input_space_sup;
% Compute input bounds for entire batch.
xBatchL = max(xBatch - epsilon,inpInf);
xBatchU = min(xBatch + epsilon,inpSup);
% Restrict noise to input bounds.
noise = 1/2*(xBatchU - xBatchL);
xBatch = 1/2*(xBatchU + xBatchL);
% Obtain generator parameters.
numInitGens = options.nn.train.num_init_gens;
initGens = options.nn.train.init_gens;
% construct input generators
if any(noise > 0)
    if strcmp(initGens,'l_inf')
        % l_inf-ball as input zonotope; the generator matrix is the 
        % indentity matrix
        xBatchG = batchG;
        idx1 = repmat(cast(1:v0,'like',xBatch),1,batchSize);
        idx2 = repelem(cast(1:batchSize,'like',xBatch),1,v0);
        % Scale non-zero entries by training perturbation noise
        nnzIdx = reshape(sub2ind(size(xBatchG),idx1,idx1,idx2), ...
            [v0 batchSize]);
        xBatchG(nnzIdx) = noise.*reshape(xBatchG(nnzIdx),[v0 batchSize]);

        if options.nn.interval_center
            xBatch = permute(cat(3,xBatch,xBatch),[1 3 2]);
        end
    elseif strcmp(initGens,'random')
        % sample a random generators from [-1,1].
        xBatchG = 2*rand(v0,numInitGens,batchSize,'like',xBatch) - 1;
        xBatchG = cast(xBatchG,'like',xBatch);
        % Scale the generators by the current training noise.
        xBatchG = epsilon*xBatchG./max(sum(abs(xBatchG),2),eps('like',xBatchG));
    elseif strcmp(initGens,'sensitivity')
        % disable backpropagation for the computation of the sensitivity
        options.nn.train.backprop = false;
        % Compute the sensitivity matrix for each input of the
        % current batch.
        [S,~] = nn.calcSensitivity(xBatch,options,false);
        % re-enable backpropagation
        options.nn.train.backprop = true;
        % Group input dimensions to have similar sensitivity.
        % Find the input pixels that affect the output the most.
        [~,idx] = sort(sum(abs(S)),2,'descend');
        % % Compute the number of entries per generator.
        % entsPerGen = floor(v0/numInitGens);
        % % Transform the indices for an efficient construction of
        % % the generator matrix.
        % idx = {':',reshape(idx(:,1:numInitGens*entsPerGen,:),1,[])};
        % % Construct the generator matrix and sum along the extra
        % % dimension to combine the unimportant generators.
        % xBatchG = reshape(idMat(idx{:}), [v0 numInitGens entsPerGen batchSize]);
        % xBatchG = reshape(sum(xBatchG,3),[v0 numInitGens batchSize]);
        xBatchG = reshape(idMat(:,idx(:,1:numInitGens,:)),[v0 numInitGens batchSize]);
        % Scale the generators by the current training noise.
        xBatchG = permute(noise,[1 3 2]).*xBatchG;

        if options.nn.interval_center
            idx_ = sub2ind([v0 batchSize],idx(:,numInitGens+1:end,:),repelem(1:batchSize',v0 - numInitGens,1));
            mask = zeros(size(xBatch),'like',xBatch);
            mask(idx_) = 1;
            xBatch = permute(cat(3,mask.*xBatchL + (1-mask).*xBatch,mask.*xBatchU + (1-mask).*xBatch),[1 3 2]);
        end
    elseif strcmp(initGens,'fgsm')
        % Compute a FGSM attack for each generator
        xsBatch = repmat(xBatch,[1 1 numInitGens]);
        xsBatchL = repmat(xBatchL,[1 1 numInitGens]);
        xsBatchU = repmat(xBatchU,[1 1 numInitGens]);
        % add random noise
        xsBatch = xsBatch + noise.*(2*rand(size(xsBatch),'like',xsBatch) - 1);
        tsBatch = repmat(tBatch,[1 1 numInitGens]);
        zBatch = nn.computePGDAttack(xsBatch(:,:),tsBatch(:,:),options, ...
            1/2*epsilon,1,xsBatchL(:,:),xsBatchU(:,:));
        zBatch = reshape(zBatch,[v0 batchSize numInitGens]);
        % Move center in the middle of attacks.
        % xBatch = 1/numInitGens*reshape(sum(zBatch,3),[v0 batchSize]);
        % Update noise.
        % noise = min(xBatchU - xBatch,xBatch - xBatchL);
        % noise = 1/2*(min(xBatch + epsilon,inpSup) - max(xBatch - epsilon,inpInf));
        % Compute generators from FGSM attacks.
        xBatchG = permute(xBatch - zBatch,[1 3 2]);
        % Scale generators.
        % s = permute(noise,[1 3 2])./sum(abs(xBatchG),2);
        % s(isnan(s)) = 0;
        % xBatchG = s.*xBatchG;
        xBatchG = 1/numInitGens*xBatchG;

        if options.nn.interval_center
            % r = reshape(sum(abs(xBatchG),2),[v0 batchSize]);
            % xBatch = interval(xBatch - (noise - r),xBatch + (noise - r));
            xBatch = permute(cat(3,xBatch,xBatch),[1 3 2]);
        end
    elseif strcmp(initGens,'patches')
        patchSize = [4 4];
        imgSize = [28 28];
        numPatch = ceil(imgSize./patchSize);

        patchImg = kron(reshape(1:prod(numPatch),numPatch),ones(patchSize,'like',xBatch));
        patchImg = patchImg(1:imgSize(1),1:imgSize(2));
        idx = sub2ind([v0 prod(numPatch)],1:v0,patchImg(:));
        patchGen = zeros([v0 prod(numPatch)],'like',xBatch);
        patchGen(idx) = 1;
        patchGen = epsilon*patchGen;

        xBatchG = repmat(patchGen,1,1,batchSize);

        if options.nn.interval_center
            xBatch = permute(cat(3,xBatch,xBatch),[1 3 2]);
        end
    else
        throw(CORAerror('CORA:wrongFieldValue', ...
            sprintf("options.nn.train.init_gens: %s!",initGens)));
    end
else
    if options.nn.interval_center
        xBatch = permute(cat(3,xBatch,xBatch),[1 3 2]);
    end
    xBatchG = zeros([v0 numInitGens batchSize],'like',xBatch);
end

% To save memory, we clear all variables that are no longer used.
batchVars = {'S','xsBatch','xsl','xsu','tsBatch','zBatch'};
clear(batchVars{:});
end

function [lossBatch] = aux_trainOtherMethods(nn,xBatch,tBatch,idxLayer,...
    options,iterScale)
% Obtain the bound of the input space.
inpInf = options.nn.train.input_space_inf;
inpSup = options.nn.train.input_space_sup;
% Handle forward and back-propagation of other training methods for robust
% neural networks, see references [3-6]
if strcmp(options.nn.train.method,'madry')
    % disable backpropagation for the computation of the pgd attacks
    options.nn.train.backprop = false;
    % Compute input bounds.
    epsilon = options.nn.train.noise;
    batchL = max(xBatch - epsilon,inpInf);
    batchU = min(xBatch + epsilon,inpSup);
    % Compute random initial start
    zInit = xBatch + epsilon*(2*rand(size(xBatch),'like',xBatch) - 1);
    % compute adversarial attack with PGD
    xBatch = nn.computePGDAttack(zInit,tBatch,options,...
        epsilon,options.nn.train.pgd_iterations,batchL,batchU, ...
        options.nn.train.pgd_stepsize,...
        options.nn.train.pgd_stepsize_decay,...
        options.nn.train.pgd_stepsize_decay_iter);
    % re-enable backpropagation
    options.nn.train.backprop = true;
    % forward propagation of perturbed inputs
    yBatch = nn.evaluate_(xBatch,options,idxLayer);
    % compute loss
    [lossBatch,grad,~] = aux_computeLoss(tBatch,yBatch,[],options);
    % backpropagation
    nn.backprop(grad,options,idxLayer);
elseif strcmp(options.nn.train.method,'gowal')     
    % linearly ramp-down \kappa: \kappa = 1,...,\kappa_{final}
    kappa = options.nn.train.kappa;
    kappa = kappa + (1 - iterScale)*(1 - kappa);
    % normal forward propagation
    yBatch = nn.evaluate_(xBatch,options,idxLayer);
    % compute regular loss
    [lossBatch,grad,~] = aux_computeLoss(tBatch,yBatch,[],options);
    % backpropagation for regular loss
    nn.backprop(kappa*grad,options,idxLayer);
    % compute input bounds
    r = options.nn.train.noise;
    % Restrict interval to input bounds.
    batchL = max(xBatch - r,inpInf);
    batchU = min(xBatch + r,inpSup);
    xBatchIval = interval(batchL,batchU);
    % elision of the last linear layer
    elision = isa(nn.layers{idxLayer(end)},'nnLinearLayer');
    if elision
        % compute output bounds (except last layer)
        yBatchIval = nn.evaluate_(xBatchIval,options,idxLayer(1:end-1));
        muk = (yBatchIval.sup + yBatchIval.inf)/2;
        rk = (yBatchIval.sup - yBatchIval.inf)/2;
        % obtain last linear layer
        linl = nn.layers{idxLayer(end)};
        % store input for backpropagation
        linl.backprop.store.input = yBatchIval;
        % compute specification matrix
        vK = size(tBatch,1);
        Ct = eye(vK,'like',tBatch) - ones(vK,1).*permute(tBatch,[3 1 2]);
        % transform weights and bias of last layer
        CtW = pagemtimes(Ct,linl.W);
        Ctb = pagemtimes(Ct,linl.b);
        % apply weights and bias of last layer
        muK = reshape(pagemtimes(CtW,permute(muk,[1 3 2])) + Ctb,size(tBatch));
        rK = reshape(pagemtimes(abs(CtW),permute(rk,[1 3 2])),size(tBatch));
        % compute robustness loss 
        zBatch = muK + rK;
        [lossBatchIval,gradIval,~] = aux_computeLoss(tBatch,zBatch,[],options);
        % backprop through last linear layer
        dW = reshape(sum(pagemtimes(Ct,'transpose',gradIval*muk','none') + ...
            pagemtimes(Ct,'transpose',gradIval*rk','none').*sign(CtW),3),size(linl.W));
        db = reshape(sum(pagemtimes(Ct,'transpose', ...
            permute(gradIval,[1 3 2]),'none'),3),size(linl.b));
        % update weights and bias
        linl.backprop.grad.('W') = linl.backprop.grad.('W') + dW;
        linl.backprop.grad.('b') = linl.backprop.grad.('b') + db;
        % backprop gradient
        dmu = reshape(pagemtimes(CtW,'transpose', ...
            permute(gradIval,[1 3 2]),'none'),size(muk));
        dr = reshape(pagemtimes(abs(CtW),'transpose', ...
            permute(gradIval,[1 3 2]),'none'),size(rk));
        % backpropagation for robustness loss
        nn.backpropIntervalBatch((1 - kappa)*(dmu - dr), ...
            (1 - kappa)*(dmu + dr),options,idxLayer(1:end-1));
    else
        % compute output bounds
        yBatchIval = nn.evaluate_(xBatchIval,options,idxLayer);
        % compute robustness loss 
        zBatch = tBatch.*yBatchIval.inf + (1-tBatch).*yBatchIval.sup;
        [lossBatchIval,gradIval,~] = aux_computeLoss(tBatch,zBatch,[],options);
        nn.backpropIntervalBatch((1 - kappa)*tBatch.*gradIval, ...
            (1 - kappa)*(1-tBatch).*gradIval,options,idxLayer);
    end
    % combine losses
    lossBatch.center = kappa*lossBatch.center + (1 - kappa)*lossBatchIval.center;
    lossBatch.vol = kappa*lossBatch.vol + (1 - kappa)*lossBatchIval.vol;
    lossBatch.total = kappa*lossBatch.total + (1 - kappa)*lossBatchIval.total;
elseif strcmp(options.nn.train.method,'trades')            
    % normal forward propagation
    yBatch = nn.evaluate_(xBatch,options,idxLayer);
    % compute regular loss
    [lossBatch,grad,~] = aux_computeLoss(tBatch,yBatch,[],options);
    % backpropagation for regular loss
    nn.backprop(grad,options,idxLayer);
    if options.nn.train.lambda > 0
        % disable backpropagation for the computation of the pgd attacks
        options.nn.train.backprop = false;
        % Compute input bounds.
        epsilon = options.nn.train.noise;
        batchL = max(xBatch - epsilon,inpInf);
        batchU = min(xBatch + epsilon,inpSup);
        % Compute random initial start
        zInit = xBatch + epsilon*(2*rand(size(xBatch),'like',xBatch) - 1);
        % compute adversarial attack with PGD
        xBatch = nn.computePGDAttack(zInit,tBatch,options,...
            epsilon,options.nn.train.pgd_iterations,batchL,batchU, ...
            options.nn.train.pgd_stepsize,...
            options.nn.train.pgd_stepsize_decay,...
            options.nn.train.pgd_stepsize_decay_iter);
        % re-enable backpropagation
        options.nn.train.backprop = true;
        % forward propagation of perturbed inputs
        zBatch = nn.evaluate_(xBatch,options,idxLayer);
        % compute robustness loss 
        [lossBatchRob,gradRob,~] = aux_computeLoss(softmax(yBatch), ...
            zBatch,[],options);
        % backpropagation for robustness loss
        nn.backprop(1/options.nn.train.lambda*gradRob,options,idxLayer);

        % combine losses
        lossBatch.center = lossBatch.center + ...
            1/options.nn.train.lambda*lossBatchRob.center;
        lossBatch.vol = lossBatch.vol + ...
            1/options.nn.train.lambda*lossBatchRob.vol;
        lossBatch.total = lossBatch.total + ...
            1/options.nn.train.lambda*lossBatchRob.total;
    end
elseif strcmp(options.nn.train.method,'sabr')     
    % clamp values to [0,1]
    batchL = max(inpInf,xBatch - options.nn.train.noise);
    batchU = min(inpSup,xBatch + options.nn.train.noise);
    % random initial value
    xBatch = batchL + (batchU - batchL).*rand(size(xBatch),'like',xBatch);
    % compute radius
    r = options.nn.train.lambda/2*(batchU - batchL);
    % determine interval center with PGD
    xBatch = nn.computePGDAttack(xBatch,tBatch,options,...
        options.nn.train.noise,options.nn.train.pgd_iterations, ...
        batchL,batchU,options.nn.train.pgd_stepsize,...
        options.nn.train.pgd_stepsize_decay,...
        options.nn.train.pgd_stepsize_decay_iter);
    % clamp center
    xBatch = min(batchU + r,max(batchL - r,xBatch));
    % compute input bounds.
    xBatchIval = interval(xBatch - r,xBatch + r);
    % compute output bounds
    yBatchIval = nn.evaluate_(xBatchIval,options,idxLayer);
    % compute robustness loss 
    zBatch = tBatch.*yBatchIval.inf + (1-tBatch).*yBatchIval.sup;
    [lossBatch,gradIval,~] = aux_computeLoss(tBatch,zBatch,[],options);
    % backpropagation for robustness loss
    nn.backpropIntervalBatch(tBatch.*gradIval,(1-tBatch).*gradIval, ...
        options,idxLayer);
else
    throw(CORAerror('CORA:wrongFieldValue', ...
        sprintf("Unsported training method: %s;",options.nn.train.method)));
end
% To save memory, we clear all variables that are no longer used.
batchVars = {'yBatch','grad','r','xBatchIval','yBatchIval','zBatch','xl','xu'};
clear(batchVars{:});
end

function aux_printParameters(options,optim,maxEpoch,miniBatchSize,maxNoise)
    % print parameters
    table = CORAtableParameters('Neural Network Training Parameters');
    table.printHeader();
    % print main parameters
    table.printContentRow('Training method',options.nn.train.method)
    table.printContentRow('Optimizer',optim.print())
    table.printContentRow('Training epochs',maxEpoch,'i')
    table.printContentRow('Mini-batch size',miniBatchSize,'i')
    table.printContentRow('Shuffle data',options.nn.train.shuffle_data)
    if strcmp(options.nn.train.method,'point')
        % no additional parameters
    else
        % set-based parameters
        table.printContentRow('(max) Training noise',maxNoise,'.2e')
        if strcmp(options.nn.train.method,'madry')
            table.printContentRow('Number of PGD iterations',options.nn.train.pgd_iterations,'i')
        elseif strcmp(options.nn.train.method,'gowal')
            table.printContentRow('kappa (gowal)',options.nn.train.kappa,'.2e')
        elseif strcmp(options.nn.train.method,'set')
            table.printContentRow('Volume heuristic for training',options.nn.train.volume_heuristic,'s')
            table.printContentRow('Weight for volume heuristic in loss (tau)',options.nn.train.tau,'.2e')
            table.printContentRow('Type of input generators',options.nn.train.init_gens,'s')
            table.printContentRow('Number of generators',options.nn.train.num_init_gens,'d')
        end
    end
    % finish table
    table.printContentRow('Early stopping after non-decr. validation steps',options.nn.train.early_stop,'i')
    table.printFooter();
end

function table = aux_setUpTrainingTable(options)
    % Specify values to print
    hvalues = {'Epoch','Iteration','Training Time',...
        sprintf('Loss (%s; train)',options.nn.train.loss),'Loss (Vol; train)',...
        'Loss (Total; train)','Loss (Total; val)'};
    % Specify formats for column values
    formats = {'d','d','s','.4e','.4e','.4e','.4e'};
    table = CORAtable('double',hvalues,formats);
end

function aux_printTrainingIteration(table,iter,epoch,trainTime,valIter,loss)

    % Set values for each column
    valIterIdx = find(valIter == iter);
    trainTimeVec = [0 0 0 0 0 trainTime];
    cvalues = {epoch,iter,datetime(trainTimeVec,'Format','HH:mm:ss'),...
        loss.center(iter),loss.vol(iter),loss.total(iter),...
        loss.valTotal(valIterIdx)};
    
    table.printContentRow(cvalues)

end

% ------------------------------ END OF CODE ------------------------------
