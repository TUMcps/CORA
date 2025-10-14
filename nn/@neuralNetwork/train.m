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
%               .propagation_batch_size: mini batch is split up to save
%                   memory
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
%        .train: total training training loss
%        .val: validation loss
%    trainTime - training time
%    
% References:
%    [1] C. M. Bishop. Pattern Recognition and Machine Learning. 
%        Information Science and Statistics. Springer New York, NY, 2006.
%    [2] Koller, L. et al. Set-based Training for Neural Network
%           Verification. TMLR 2025
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
%                16-January-2024 (sensitivity-based generators)
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
    {trainX,'att','numeric'}; ...
    {trainT,'att','numeric'}; ... 
    {valX,'att','numeric'}; ... 
    {valT,'att','numeric'}; ... 
    {options,'att','struct'}; ... 
    {verbose,'att','logical'};
});
% Set default training parameters
options = nnHelper.validateNNoptions(options,true);

% Extract training parameters
maxEpoch = options.nn.train.max_epoch;
miniBatchSize = options.nn.train.mini_batch_size;
propBatchSize = min(miniBatchSize,options.nn.train.propagation_batch_size);

maxNoise = options.nn.train.noise;
numWarmUpEpochs = options.nn.train.warm_up;
% Compute numbre of ramp up epochs
numRampUpEpochs = options.nn.train.ramp_up - numWarmUpEpochs;

% Obtain number of input dimension and number of data points
[~,N] = size(trainX);
[v0,~] = size(valX);
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
    'train',zeros(numIter,1),...
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
end

% set up verbose output table
table = aux_setUpTrainingTable(options);
if verbose
    table.printHeader();
end

% Start training time
timerVal = tic;

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
        iterScale = (iter - numWarmUpEpochs*iterPerEpoch)...
            /(numRampUpEpochs*iterPerEpoch);
        iterScale = min(max(0,iterScale),1);
        loss.iterScale(iter) = iterScale;
        % Scale training noise
        options.nn.train.noise = maxNoise*iterScale; % *noiseScale;

        % Compute indices for elements in current batch
        miniBatchIdx = i:i+min(miniBatchSize,N)-1;
        % Extract data of current batch
        xBatch = cast(trainX(:,miniBatchIdx),'like',inputDataClass);
        tBatch = cast(trainT(:,miniBatchIdx),'like',inputDataClass);
        % Data augmentation.
        if isfield(options.nn.train,'data_augmentation')
            xBatch = options.nn.train.data_augmentation(xBatch);
        end
        % enable backpropagation
        options.nn.train.backprop = true;
        options.nn.interval_center = ivalC;
        if strcmp(options.nn.train.method,'point') ...
                || options.nn.train.noise == 0 ... (strcmp(options.nn.train.method,'set') && options.nn.train.noise == 0)
            % Forward propagation
            yBatch = nn.evaluate_(xBatch,options,idxLayer);
            % Compute loss
            [lossBatch,grad,~] = aux_computeLoss(tBatch,yBatch,[],options);
            % Backpropagation
            nn.backprop(grad,options,idxLayer);
        elseif strcmp(options.nn.train.method,'set')
            % Split mini-batch and propagate parts to save memory.
            lossBatch.center = 0;
            lossBatch.vol = 0;
            lossBatch.train = 0;
            % lossBatch.approxErrRel = 0;
            lossBatch.minMeanMaxVarGG = 0;
            lossBatch.minMeanMaxVarGR = 0;
            lossBatch.minMeanMaxVarGc = 0;
            for j=1:propBatchSize:...
                    max(1,miniBatchSize - mod(miniBatchSize,propBatchSize))
                % Compute indices for elements in current batch
                propIdx = j:j+min(propBatchSize,miniBatchSize)-1;

                % Batch norm is based on statistics of the nominal input; normal forward
                % propagation to compute stats.
                options.nn.batch_norm_calc_stats = true; % compute batch stats.

                % construct input generators
                [xBatch_,GBatch_] = aux_construtInputGenerators(nn, ...
                    xBatch(:,propIdx),... idBatch(:,:,1:length(propIdx)), ...
                        idMat,tBatch(:,propIdx),options);
                if strcmp(options.nn.approx_error_order,'sensitivity*length')
                    % Compute sensitivity; need to select most influential
                    % approximation errors.
                    nn.calcSensitivity(xBatch(:,propIdx));
                end
                % Set-based forward propagation.
                [cyBatch,GyBatch] = nn.evaluateZonotopeBatch_( ...
                    xBatch_,GBatch_,options,idxLayer);
                % Compute set-based loss.
                [lossBatch_,gc,gG] = aux_computeLoss(tBatch(:,propIdx), ...
                    cyBatch,GyBatch,options);
                % Set-based backpropagation.
                nn.backpropZonotopeBatch_(gc,gG,options,idxLayer,true);
                % Add propgation loss.
                lossScale = length(propIdx)/length(miniBatchIdx);
                lossBatch.center = lossBatch.center + ...
                    lossScale*lossBatch_.center;
                lossBatch.vol = lossBatch.vol + ...
                    lossScale*lossBatch_.vol;
                lossBatch.train = lossBatch.train + ...
                    lossScale*lossBatch_.train;
            end
        else
            % Train using other methods ('madry','gowal','trades','sabr'), 
            % see [3-6]
            lossBatch = aux_trainOtherMethods(nn,xBatch,tBatch,idxLayer,...
                options,iterScale);
        end
        % Store loss of current batch
        loss.center(iter) = loss.center(iter) + lossBatch.center;
        loss.vol(iter) = loss.vol(iter) + lossBatch.vol;
        loss.train(iter) = loss.train(iter) + lossBatch.train;
        % loss.noiseScale(iter) = noiseScale;

        % To save memory, we clear all variables that are only used
        % during the batch computation
        batchVars = {'xBatch','tBatch','yBatch','grad','GBatch','cyBatch',...
            'GyBatch','gc','gG'};
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
            % Loop index.
            k = 1;
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
                lossValCenter(k) = lossValBatch.center*length(valBatchIdx);
                lossValVol(k) = lossValBatch.vol*length(valBatchIdx);
                lossValTotal(k) = lossValBatch.train*length(valBatchIdx);
                % Increment loop index.
                k = k + 1;

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
            trainTime = toc(timerVal);
            aux_printTrainingIteration(options,table,iter,epoch,trainTime,valIter,loss)
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
    % Check if we stopped early
    if numNonDecrVal > options.nn.train.early_stop
        break;
    end
end

% print bottom
if verbose
    table.printFooter();
end

% Clear variables that are no longer used.
vars = {'idMat'}; ... 'idBatch'
clear(vars{:});

% (potentially) gather weights of the network from gpu
nn.castWeights(single(1));

% Delete gradients
optim.deleteGrad(nn,options);

end


% Auxiliary functions -----------------------------------------------------

function [loss,gc,gG] = aux_computeLoss(t,cy,Gy,options)
% Compute the loss for a given zonotope Y=(cy,Gy). If the generator matrix
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
    % cy = cy - max(cy,[],1);
    subMax = @(x) x - max(x,[],1);
    % Use numerically stable log-softmax
    logSoftmax = @(x) x - log(sum(exp(x),1));
    % Classification loss (cross-entropy)
    lossFun = @(t,y) -1/batchSize*sum(t.*logSoftmax(subMax(y)),'all');
    gradFun = @(t,y) 1/batchSize*(softmax(subMax(y)) + (-t));
elseif strcmp(options.nn.train.loss,'custom')
    % Take given loss function and derivative
    lossFun = options.nn.train.loss_fun;
    gradFun = options.nn.train.loss_grad;
    % Compute custom loss
    loss.train = lossFun(t,cy,Gy);
    loss.center = 0;
    loss.vol = 0;
    [gc,gG] = gradFun(t,cy,Gy);
    return
else
    throw(CORAerror('CORA:wrongFieldValue',...
        'options.nn.train.loss',{'mse','softmax+log','custom'}));
end

if options.nn.interval_center && ndims(cy) == 3
    cyl = reshape(cy(:,1,:),[vK batchSize]);
    cyu = reshape(cy(:,2,:),[vK batchSize]);

    % Center Loss
    z = 1/2*(cyu + cyl);
    % Compute loss of center
    loss.center = cast(lossFun(t,z),'like',single(1));  
   
    % Compute gradient of center
    gc = gradFun(t,z);
    gc = 1/2*permute(cat(3,gc,gc),[1 3 2]);
   
    idMat = eye(size(z,1),'like',z);
    Gy = cat(2,Gy,1/2*permute(cyu - cyl,[1 3 2]).*idMat);
else
    % Use center.
    z = cy;
    % Compute loss of center.
    loss.center = cast(lossFun(t,z),'like',single(1));
    % Compute gradient of center.
    gc = gradFun(t,z);
end


if isempty(Gy) 
    % Training point-based: There is no volume loss
    loss.vol = 0;
    loss.train = loss.center;
    % There is not gradient for the generator matrix
    gG = zeros(size(Gy),'like',Gy);
else 
    % Compute volume heuristic
    if strcmp(options.nn.train.volume_heuristic,'interval')
        % The interval norm is the sum of all absolute values of the
        % generator matrix.
        yPredVolHeu = 1/batchSize*sum(abs(Gy),'all');
        gG = 1/batchSize*sign(Gy);
    elseif strcmp(options.nn.train.volume_heuristic,'f-radius')
        % The F-radius is the square-root of the sum of all squared entries
        % of the generator matrix.
        frad = sqrt(sum(Gy.^2,[1 2]));
        zIdx = (frad(:) < eps('like',frad));
        yPredVolHeu = 1/batchSize*sum(frad);
        gG = 1/batchSize*(Gy./frad);
        gG(:,:,zIdx) = 0;
    elseif strcmp(options.nn.train.volume_heuristic,'set-eval')
        % Set-based evaluation of the gradient.
        yPredVolHeu = 1/batchSize*1/2*sum(Gy.^2,'all');
        gG = 1/batchSize*Gy;
    else
        throw(CORAerror('CORA:wrongFieldValue',...
            'options.nn.train.volume_heuristic',{'interval','f-radius'}));
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
    loss.train = cCenter*loss.center + cVol*loss.vol;
    % Scale the gradients according to their weights
    gc = cCenter*gc;
    gG = cVol*gG; 

    if options.nn.interval_center % && ~options.nn.train.ibp_center_loss
        % Extract gradient for the interval center.
        q = size(Gy,2) - vK;
        diagIdx = reshape(sub2ind(size(Gy),repmat(1:vK,1,batchSize),...
            q + repmat(1:vK,1,batchSize),repelem(1:batchSize,1,vK)),...
            [vK batchSize]);
        gc = gc + 1/2*permute(cat(3,-gG(diagIdx),gG(diagIdx)),[1 3 2]);
        gG = gG(:,1:q,:);
    end
end

end

function [xBatch,GBatch] = ...
    aux_construtInputGenerators(nn,xBatch,idMat,tBatch,options) ... ,idBatch
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
        GBatch = permute(noise,[1 3 2]).*idMat;

        if options.nn.interval_center
            xBatch = permute(cat(3,xBatch,xBatch),[1 3 2]);
        end
    elseif strcmp(initGens,'random')
        % sample a random generators from [-1,1].
        GBatch = 2*rand(v0,numInitGens,batchSize,'like',xBatch) - 1;
        GBatch = cast(GBatch,'like',xBatch);
        % Scale the generators by the current training noise.
        GBatch = epsilon*GBatch./max(sum(abs(GBatch),2),eps('like',GBatch));
    elseif strcmp(initGens,'sensitivity')
        
        % 1. Forward propagation, with backprop enabled.
        options.nn.train.backprop = true;
        yBatch = nn.evaluate(xBatch,options);
        % 2. Backpropagation.
        idxLayer = 1:length(nn.layers);
        gradBatch = nn.backprop(yBatch - tBatch,options,idxLayer,false);
        % Find the input pixels that affect the output the most.
        [~,dimIdx] = sort(gradBatch,'descend');

        % Initialize generator matrix.
        GBatch = zeros([v0 numInitGens batchSize],'like',xBatch);
        % Compute indices for non-zero entries.
        gIdx = sub2ind(size(GBatch),...
            dimIdx(1:numInitGens,:),...
            repmat((1:numInitGens)',1,batchSize),...
            repelem(1:batchSize,numInitGens,1));
        % Permute noise vector.
        noise_ = noise(dimIdx);
        % Set non-zero generator entries.
        GBatch(gIdx) = noise_(1:numInitGens,:);

        if options.nn.interval_center
            % Move remaining noise to the interval center.
            rBatch = sum(GBatch,2);
            xBatch = permute(cat(3,xBatchL + rBatch,xBatchU - rBatch),[1 3 2]);
        end
    elseif startsWith(initGens,'fgsm')
        % Compute a FGSM attack for each generator
        xsBatch = repmat(xBatch,[1 1 numInitGens]);
        xsBatchL = repmat(xBatchL,[1 1 numInitGens]);
        xsBatchU = repmat(xBatchU,[1 1 numInitGens]);        
        % add random noise
        xsBatch = xsBatch + noise.*(2*rand(size(xsBatch),'like',xsBatch) - 1);
        tsBatch = repmat(tBatch,[1 1 numInitGens]);
        zBatch = nn.computePGDAttack(xsBatch(:,:),tsBatch(:,:),options, ...
            epsilon,1,xsBatchL(:,:),xsBatchU(:,:), ...
            epsilon,0,[]); % @(t,y) 1 - t);
        zBatch = reshape(zBatch,[v0 batchSize numInitGens]);

        % Move center in the middle of attacks.
        xBatch = 1/numInitGens*reshape(sum(zBatch,3),[v0 batchSize]);
        % xBatch = 1/2*reshape(max(zBatch,[],3) + min(zBatch,[],3),[v0 batchSize]);
        % Update noise.
        noise = min(xBatchU - xBatch,xBatch - xBatchL);

        % Scale generators.
        if endsWith(initGens,'1')
            % Compute generators from FGSM attacks.
            GBatch = permute(xBatch - zBatch,[1 3 2]);
            s = permute(noise,[1 3 2])./sum(abs(GBatch),2);
            s(isnan(s)) = 0;
            GBatch = s.*GBatch;
        elseif endsWith(initGens,'2')
            % Compute generators from FGSM attacks.
            GBatch = 1/numInitGens*permute(xBatch - zBatch,[1 3 2]);
        elseif endsWith(initGens,'3')
            GBatch = zeros([v0 numInitGens batchSize],'like',xBatch);
        end

        if options.nn.interval_center
            r = reshape(sum(abs(GBatch),2),[v0 batchSize]);
            xBatch = permute(cat(3,xBatch - (noise - r),xBatch + (noise - r)),[1 3 2]);
            % xBatch = permute(cat(3,xBatch,xBatch),[1 3 2]);
        end
    else
        throw(CORAerror('CORA:wrongFieldValue', ...
            'options.nn.train.init_gens',{'l_inf','random','sensitivity','fgsm'}));
    end
else
    if options.nn.interval_center
        xBatch = permute(cat(3,xBatch,xBatch),[1 3 2]);
    end
    GBatch = zeros([v0 numInitGens batchSize],'like',xBatch);
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

% Batch norm is based on statistics of the nominal input; normal forward
% propagation to compute stats.
options.nn.batch_norm_calc_stats = true; % compute batch stats.
options.nn.train.backprop = false; % Disable backprop.
nn.evaluate_(xBatch,options,idxLayer);
% Reset flags.
options.nn.batch_norm_calc_stats = false;
options.nn.train.backprop = true;
options.nn.batch_norm_moving_stats = true;

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
    elision = false; % isa(nn.layers{idxLayer(end)},'nnLinearLayer');
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
    lossBatch.train = kappa*lossBatch.train + (1 - kappa)*lossBatchIval.train;
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
        lossBatch.train = lossBatch.train + ...
            1/options.nn.train.lambda*lossBatchRob.train;
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
    throw(CORAerror('CORA:wrongFieldValue',...
        'options.nn.train.method',{'point','set','mardy','gowal','trades','sabr'}));
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
    % Specify values to print and formats for column values
    if strcmp(options.nn.train.method,'set') && ~strcmp(options.nn.train.loss,'custom')
        hvalues = {'Epoch','Iteration','Training Time',...
            sprintf('Loss (%s; train)',options.nn.train.loss),...
                'Loss (Vol; train)','Loss (Total; train)','Loss (Total; val)'};
        formats = {'d','d','s','.4e','.4e','.4e','.4e'};
    else
        hvalues = {'Epoch','Iteration','Training Time',...
            'Loss (train)','Loss (val)'};
        formats = {'d','d','s','.4e','.4e'};
    end
    table = CORAtable('double',hvalues,formats);
end

function aux_printTrainingIteration(options,table,iter,epoch,trainTime,valIter,loss)

    % Set values for each column
    valIterIdx = find(valIter == iter);
    trainTimeVec = [0 0 0 0 0 trainTime];
    if strcmp(options.nn.train.method,'set') && ~strcmp(options.nn.train.loss,'custom')
        cvalues = {epoch,iter,datetime(trainTimeVec,'Format','HH:mm:ss'),...
            loss.center(iter),loss.vol(iter),loss.train(iter),...
                loss.valTotal(valIterIdx)};
    else
        cvalues = {epoch,iter,datetime(trainTimeVec,'Format','HH:mm:ss'),...
            loss.train(iter),loss.valTotal(valIterIdx)};
    end
    
    table.printContentRow(cvalues)

end

% ------------------------------ END OF CODE ------------------------------
