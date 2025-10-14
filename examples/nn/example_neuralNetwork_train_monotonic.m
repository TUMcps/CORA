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
%    [1] Kitouni, O. et al. Expressive monotonic neural networks. (ICLR). 2023
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

% We train a monotonic neural network to approximate f(x) = pi*x + sin(pi*x).

f = @(x) x + 1/(2*pi)*sin(2*pi*x);

% Generate the Dataset. ---------------------------------------------------

% Specify the bounds of the input space.
xl = 0;
xu = 1;

% Samples data points.
xs = linspace(xl,xu,1000);
ts = f(xs);

% Specify Network and Training Parameters. --------------------------------

n0 = 1; % number of input dimensions.
nK = 1; % number of output dimensions.
nk = 50; % number of hidden dimensions.
K = 3; % number of layers.
actfun = 'groupSort'; % type of activation function, i.e., {'relu','tanh','groupSort'}.

% Specify the training parameters.
lr = 1e-3; % learning rate.
numEpoch = 300; % number of epochs.
bSz = 64; % batch size.

% Train the Neural Network. -----------------------------------------------

% Specify training options.
options.nn.train = struct( ...
    'optim',nnAdamOptimizer(lr),...
    'max_epoch',numEpoch,...
    'mini_batch_size',bSz,...
    'loss','mse',...
    'shuffle_data','every_epoch' ...
);

% Create a random neural network.
lambda = 1; % Lipschitz constant
nn = aux_generateMonotonicNeuralNetwork(options,n0,nk,nK,K,actfun,lambda);

% Train the neural network.
loss = nn.train(xs,ts,[],[],options,true);

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

function [nn,layers] = aux_generateMonotonicNeuralNetwork(options, ...
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
        % Create a linear layer.
        layers{end+1} = nnLipConstrLinearLayer(zeros(nout,nin),zeros(nout,1),lambda^-(1/K));
        if i < length(nks)-2
            % Append an activation layer (only if not the last linear layer).
            switch actfun
                case 'relu'
                    layers{end+1} = nnReLULayer;
                case 'tanh'
                    layers{end+1} = nnTanhLayer;
                case 'groupSort'
                    layers{end+1} = nnGroupSortLayer;
                otherwise
                    throw(CORAerror('CORA:wrongValue', ...
                        'actfun',{'relu','tanh'}));
            end
        end
    end
    % Create the resiudal connection.
    addLayer = nnCompositeLayer({ ...
        {nnLinearLayer(lambda*ones(n0,1),0,'residualConnection',false)}; ... % Linear layer to sum input dimensions with fixed weights.
        layers ...
    },'add');

    % Initialize network.
    nn = neuralNetwork({addLayer});
    nn.setInputSize([n0 1]);
    % Initialize the weights and bias.
    nn.initWeights('glorot');
    % Normalize the weights.
    nn.normWeights(options);

end

% ------------------------------ END OF CODE ------------------------------
