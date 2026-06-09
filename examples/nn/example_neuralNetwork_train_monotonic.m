function res = example_neuralNetwork_train_monotonic()
% example_neuralNetwork_train_monotonic - example for training a monotonic 
%   neural network.
%
% Syntax:
%    res = example_neuralNetwork_train_monotonic()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% References:
%    [1] Nolte, N. et al. Expressive monotonic neural networks. (ICLR). 2023
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       18-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng('default')

% We train a monotonic neural network to approximate h.
a = 3/2;
f = @(x) x + sin(x); % Base monotonic function
g = @(x) f(a*2*pi*x); % Scale the input to [0,1]
h = @(x) 1/g(1)*g(x); % Scale the output to [0,1]

% Generate the Dataset. ---------------------------------------------------

% Specify the bounds of the input space.
xl = 0;
xu = 1;

% Samples data points.
xs = linspace(xl,xu,1000);
ts = h(xs);
% Create validation data.
vXs = linspace(xl,xu,333);
vTs = h(vXs);

% Specify Network and Training Parameters. --------------------------------

n0 = size(xs,1); % number of input dimensions.
nK = size(ts,1); % number of output dimensions.
nk = 32; % number of hidden dimensions.
K = 4; % number of layers.
actfun = 'groupsort'; % type of activation function

% Specify the training parameters.
lr = 1e-2; % learning rate.
numEpoch = 250; % number of epochs.
bSz = 64; % batch size.

% Train the Neural Network. -----------------------------------------------

% Specify training options.
options.nn.train = struct( ...
    'optim',nnAdamOptimizer(lr),...
    'max_epoch',numEpoch,...
    'mini_batch_size',bSz,...
    'loss','mse' ...
);

% Create a random neural network.
lambda = 3; % Lipschitz constant
[nn,~,options] = aux_generateMonotonicNeuralNetwork(options,n0,nk,nK,K,actfun,lambda);

% Train the neural network.
loss = nn.train(xs,ts,vXs,vTs,options,true);

% Visualize the Loss. -----------------------------------------------------

figure; hold on;
title('Training Loss')
xlabel('#Epoch')
ylabel('Training Loss')
plot(1:length(loss.train),loss.train,'DisplayName','Training Loss');
legend

% Visualize the trained Neural Network. -----------------------------------

% Compute the output of the neural network.
ys = nn.evaluate(xs);

figure; hold on;
title('Neural Network Approximation')
xlabel('Input')
ylabel('Output')
plot(xs,ts,'DisplayName','Target Function');
plot(xs,ys,'DisplayName','Neural Network');
legend

res = true;

end


% Auxiliary functions -----------------------------------------------------

function [nn,layers,options] = aux_generateMonotonicNeuralNetwork(options, ...
    n0,nk,nK,K,actfun,lambda)
    % Generate a random neural network.
    % - n0: number of input dimensions.
    % - nk: number of hidden dimensions.
    % - nK: number of output dimensions.
    % - K: number of layers.
    % - actfun: type of activation function, 
    %       i.e., {'relu','tanh','groupSort'}.
    % - lambda: Lipschitz constant

    % Set default options parameters.
    options = nnHelper.validateNNoptions(options,true);
    
    % Compute number of neurons per linear layer.
    nks = [n0 repmat(nk,1,K) nK];

    % Initialize a cell array to store the layers.
    layers = {};
    % Create the layers.
    for i=1:length(nks)-1
        % Obtain the number of input neurons.
        nin = nks(i);
        % Obtain the number of output neurons.
        nout = nks(i+1);
        % Compute the Lipschitz constant for the i-th layer.
        lambdai = nthroot(lambda,K);
        % Create a linear layer.
        layers{end+1} = nnLipConstrLinearLayer( ...
            zeros(nout,nin),zeros(nout,1),lambdai);
        if i < length(nks)-2
            % Create an activation layer.
            actl = nnActivationLayer.instantiateFromString(actfun);
            % Append the activation layer.
            layers{end+1} = actl;
        end
    end
    % Create the resiudal connection.
    addLayer = nnCompositeLayer({ ...
        {nnLinearLayer(lambda*ones([nK n0]),0,'residualConnection',false)}; ... % Linear layer to sum input dimensions with fixed weights.
        layers ...
    },'add');

    % Initialize network.
    nn = neuralNetwork({addLayer});
    nn.setInputSize([n0 1]);
    % Initialize the weights and bias.
    nn.initWeights('glorot');
end

% ------------------------------ END OF CODE ------------------------------
