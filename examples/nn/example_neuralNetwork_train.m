function res = example_neuralNetwork_train()
% example_neuralNetwork_train - example for training a neural network.
%
% Syntax:
%    res = example_neuralNetwork_train()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       18-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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
maxEpoch = 200;
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

% Initialize network.
nn = neuralNetwork(layers);
nn.initWeights();

% Train Networks ----------------------------------------------------------

options.nn = struct(...
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
loss = nn.train(trainX,trainY,valX,valY,options,true);

% Compute Output Bounds ---------------------------------------------------

Y = nn.evaluate(zonotope(1/2*(xu + xl), 1/2*diag(xu - xl))); 
l = Y.c - sum(abs(Y.G),2);
u = Y.c + sum(abs(Y.G),2);

res = true;

end

% ------------------------------ END OF CODE ------------------------------
