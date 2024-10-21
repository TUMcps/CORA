function res = testnn_neuralNetwork_train()
% testnn_neuralNetwork_train - unit test function for 
%     neuralNetwork/train: train a neural network with CORA and the
%     DeepLearning Toolbox with same training parameters and compare the 
%     weights (Regression Task).
%
% Syntax:
%    res = testnn_neuralNetwork_train()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/train

% Authors:       Lukas Koller
% Written:       05-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Load Dataset ------------------------------------------------------------

% Specify input bounds.
xl = -pi;
xu = pi;

% learn to draw a circle
trainX = linspace(xl,xu,500);
trainY = [sin(trainX') cos(trainX')]';

valX = linspace(xl,xu,133);
valY = [sin(valX') cos(valX')]';

% Set Network and Training Parameters -------------------------------------

% network parameters
numHiddenLayers = 3;
numHiddenNeurons = 40;

% retrieve input and output dimensions
dimIn = size(trainX,1);
dimOut = size(trainY,1);

% compute number of neurons per linear layer
numNeurons = [dimIn repmat(numHiddenNeurons,1,numHiddenLayers) dimOut];

% training parameters
lr = 1e-2; % learning rate
maxEpoch = 100;
beta = 0.9; % momentum
miniBatchSize = 64;

% Initialize Network Parameters -------------------------------------------

% construct layers[layers; nnTanhLayer]
dltLayers = [featureInputLayer(dimIn)];
for i=2:length(numNeurons)-1
    outputSize = numNeurons(i);
    fcLayer = fullyConnectedLayer(outputSize);
    dltLayers = [dltLayers; fcLayer; tanhLayer];
end
dltLayers = [dltLayers; fullyConnectedLayer(numNeurons(end)); sigmoidLayer];

% initialize weights and biases of the layers
nnDltInit = dlnetwork(dltLayers);

% Train Network (DEEP-LEARNING TOOLBOX) -----------------------------------

% set training options
dltOptions = trainingOptions('sgdm',... % stochastic Gradient-Descent optimizer
    MaxEpochs=maxEpoch,...
    Momentum=beta,...
    InitialLearnRate=lr,...
    MiniBatchSize=miniBatchSize,...
    Shuffle='never',...
    L2Regularization=0.0,...
    Verbose=false,...
    ValidationData={valX',valY'}...
);
% for trainNetwork we need to append an output layer
dltLayers = [nnDltInit.Layers; regressionLayer];
% train neural network
nnDlt = trainNetwork(trainX',trainY',dltLayers,dltOptions);

% Train Network (CORA) ----------------------------------------------------

% convert initialized but un-trained DL-toolbox NN to CORA-NN
layers = num2cell(nnDltInit.Layers);
nnCora = neuralNetwork.convertDLToolboxNetwork(layers,false);

% instantiate training options
options.nn = struct(...
    'poly_method','bounds',...
    'train',struct( ...
        'optim',nnSGDOptimizer(lr,beta),...
        'max_epoch',maxEpoch,...
        'mini_batch_size',miniBatchSize,...
        'loss','mse',...
        'method','point'...
    )...
);

loss = nnCora.train(trainX,trainY,valX,valY,options,false);

% Compare the Network Weights ---------------------------------------------

% accumulate the difference between the weights and biases
weightsDiff = [];
for i=1:length(nnCora.layers)
    coraLayeri = nnCora.layers{i};
    if isa(coraLayeri,'nnLinearLayer')
        dltLayeri = nnDlt.Layers(i+1);
        weightsDiffi = [coraLayeri.W - dltLayeri.Weights,coraLayeri.b - dltLayeri.Bias];
        weightsDiff = [weightsDiff reshape(weightsDiffi,1,[])];
    end
end

assert(all(abs(weightsDiff) <= 1e-6)); % Deep-Learning Toolbox uses single

% Check Edge Cases for Training Parameters --------------------------------

% Batch size larger than training data.
loss = nnCora.train(trainX(:,1:miniBatchSize-1),...
    trainY(:,1:miniBatchSize-1),valX(:,1:miniBatchSize-1),...
    valY(:,1:miniBatchSize-1),options,true);
% One iteration per epoch.
numIterations = maxEpoch;

assert(length(loss.center) == numIterations);

% test completedte
res = true;

% ------------------------------ END OF CODE ------------------------------
