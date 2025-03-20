function res = testnn_neuralNetwork_train_set()
% testnn_neuralNetwork_train_set - unit test function for 
%     neuralNetwork/train: train a neural network set-based with CORA (with
%     epsilon=0); compare with point-based training.
%
% Syntax:
%    res = testnn_neuralNetwork_train_set()
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
% Written:       07-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

rng('default')

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

layers = {};
for i=1:length(numNeurons)-2
    nin = numNeurons(i);
    nout = numNeurons(i+1);
    layers{end+1} = nnLinearLayer(zeros(nout,nin),zeros(nout,1));
    layers{end+1} = nnTanhLayer;
end
nin = numNeurons(end-1);
nout = numNeurons(end);
layers{end+1} = nnLinearLayer(zeros(nout,nin),zeros(nout,1));

% Initialize point-based network.
nnP = neuralNetwork(layers);
nnP.initWeights();

% Copy the point-based network.
nnS = nnP.copyNeuralNetwork();

% Train Networks ----------------------------------------------------------

% instantiate training options: point and set-based with not perturbation
optionsP.nn = struct(...
    'poly_method','bounds',...
    'train',struct( ...
        'optim',nnSGDOptimizer(lr,beta),...
        'max_epoch',maxEpoch,...
        'mini_batch_size',miniBatchSize,...
        'loss','mse',...
        'method','point'...
    )...
);
% Train point-based.
lossP = nnP.train(trainX,trainY,valX,valY,optionsP,false);

% Specify the training options.
optionsS.nn = struct(...
    'use_approx_error',true,...
    'poly_method','bounds',...
    'train',struct( ...
        'optim',nnSGDOptimizer(lr,beta),...
        'max_epoch',maxEpoch,...
        'mini_batch_size',miniBatchSize,...
        'loss','mse',...
        'method','set',...
        'noise',0.5,...
        'input_space_inf',xl,...
        'input_space_sup',xu,...
        'tau',0.5,...
        'volume_heuristic','f-radius',...
        'zonotope_weight_update','sum',...
        'num_approx_err',inf,...
        'init_gens','l_inf',...
        'num_init_gens',inf...
    )...
);
% Train set-based.
lossS = nnS.train(trainX,trainY,valX,valY,optionsS,false);

% Compare Output Bounds ---------------------------------------------------

Yp = nnP.evaluate(zonotope(1/2*(xu + xl), 1/2*diag(xu - xl))); 
lp = Yp.c - sum(abs(Yp.G),2);
up = Yp.c + sum(abs(Yp.G),2);

Ys = nnS.evaluate(zonotope(1/2*(xu + xl), 1/2*diag(xu - xl))); 
ls = Ys.c - sum(abs(Ys.G),2);
us = Ys.c + sum(abs(Ys.G),2);

% Verify that set-based bounds are tighter.
assert(all(lp < ls))
assert(all(us < up));

% ------------------------------ END OF CODE ------------------------------
