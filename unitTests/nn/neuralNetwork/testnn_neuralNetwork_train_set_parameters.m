function res = testnn_neuralNetwork_train_set_parameters()
% testnn_neuralNetwork_train_set_parameters - unit test function for 
%     neuralNetwork/train: check different parameters
%
% Syntax:
%    res = testnn_neuralNetwork_train_set_parameters()
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
% Written:       08-April-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

rng('default')

% Load Dataset ------------------------------------------------------------

% Specify input bounds.
xl = 0;
xu = 1;

% learn to draw a circle
trainX = linspace(0,1,1000);
trainY = [sin(2*pi*(trainX' - 1)) cos(2*pi*(trainX' - 1))]';

valX = linspace(xl,xu,133);
valY = [sin(2*pi*(valX' - 1)) cos(2*pi*(valX' - 1))]';

% Set Network and Training Parameters -------------------------------------

% network parameters
numHiddenLayers = 3;
numHiddenNeurons = 50;

% retrieve input and output dimensions
dimIn = size(trainX,1);
dimOut = size(trainY,1);

% compute number of neurons per linear layer
numNeurons = [dimIn repmat(numHiddenNeurons,1,numHiddenLayers) dimOut];

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

% Initialize neural network network.
nn = neuralNetwork(layers);
nn.initWeights();

% Train Networks ----------------------------------------------------------

% Specify the training options.
options.nn = struct(...
    'use_approx_error',true,...
    'poly_method','bounds',...
    'train',struct( ...
        'optim',nnAdamOptimizer,...
        'max_epoch',200,...
        'mini_batch_size',128,...
        'shuffle_data','every_epoch',...
        'loss','mse',...
        'method','set',...
        'noise',1e-3,...
        'tau',1e-3,...
        'input_space_inf',xl,...
        'input_space_sup',xu,...
        'volume_heuristic','f-radius',...
        'zonotope_weight_update','sum',...
        'exact_backprop',true...
    )...
);
options.nn.train.init_gens = 'l_inf';
options.nn.train.num_init_gens = inf;
options.nn.train.num_approx_err = inf;
options.nn.train.approx_error_order = 'sensitivity*length';
% Train set-based.
loss = nn.train(trainX,trainY,valX,valY,options,true);
% Check the output of the neural network.
assert(all(withinTol(nn.evaluate(valX),valY,3e-1),'all'));

% Reset weights.
nn.initWeights();
% Modify options and retrain.
options.nn.train.warm_up = 10;
options.nn.train.ramp_up = 50;
options.nn.train.init_gens = 'fgsm2';
options.nn.train.num_init_gens = 5;
options.nn.train.num_approx_err = 49;
options.nn.train.approx_error_order = 'length';
% Train set-based.
loss = nn.train(trainX,trainY,valX,valY,options,true);
% Check the output of the neural network.
assert(all(withinTol(nn.evaluate(valX),valY,3e-1),'all'));

% Reset weights.
nn.initWeights();
% Modify options and retrain.
options.nn.interval_center = true;
% Train set-based.
loss = nn.train(trainX,trainY,valX,valY,options,true);
% Check the output of the neural network.
assert(all(withinTol(nn.evaluate(valX),valY,3e-1),'all'));

% Reset weights.
nn.initWeights();
% Modify options and retrain.
options.nn.train.use_gpu = false;
options.nn.train.exact_backprop = false;
% Train set-based.
loss = nn.train(trainX,trainY,valX,valY,options,true);
% Check the output of the neural network.
assert(all(withinTol(nn.evaluate(valX),valY,3e-1),'all'));

% ------------------------------ END OF CODE ------------------------------
