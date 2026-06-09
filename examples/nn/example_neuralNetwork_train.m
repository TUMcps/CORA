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

% Generate two moons dataset.
N = 1000; % number of samples
[xs,ts] = aux_generateTwoMoons(N,0.25);
% Specify input bounds.
xl = [-1.5; -1.5];
xu = [2.5; 1.5];
% Generate validation data.
[vXs,vTs] = aux_generateTwoMoons(floor(N/3),0.25);

% Set Network and Training Parameters -------------------------------------

% Specify the network parameters.
n0 = size(xs,1); % number of input dimensions.
nK = size(ts,1); % number of output dimensions.
nk = 16; % number of neurons per hidden layer.
K = 2; % number of hidden linear layers.
actfun = 'tanh'; % type of activation function, e.g., 'relu','tanh',...

% Specify the main training parameters.
lr = 5e-3; % learning rate
maxEpoch = 50; % total number of training epochs
bSz = 64; % batch size

% Train a Neural Network --------------------------------------------------

% Specify training options.
options.nn = struct(...
    'use_approx_error',true,...
    'poly_method','bounds',...
    'train',struct( ...
        'optim',nnAdamOptimizer(lr),...
        'max_epoch',maxEpoch,...
        'mini_batch_size',bSz,...
        'loss','softmax+log',...
        ... Specify the parameter for set-based training
        'method','set',... training method
        'noise',0.25, ... training noise
        'input_space_inf',xl,... lower bound of the input space 
        'input_space_sup',xu,... upper bound of the input space
        'tau',0.1 ... set-based loss weighting parameter
    )...
);

% Initialize the Neural Network. ------------------------------------------
[nn,~,options] = aux_generateNeuralNetwork(options,n0,nk,nK,K,actfun);

% Train set-based.
loss = nn.train(xs,ts,vXs,vTs,options,true);

% Compute the accurary on the validation data. Compute output samples.
vYs = nn.evaluate(vXs);
% Obtain class labels.
[~,vl] = max(vYs,[],1);
% Obtain predicted classes.
[~,vk] = max(vYs,[],1);
fprintf('Accuracy: %.2f\n',sum(vk == vl)/length(vl)*100);

% Visualize the Loss. -----------------------------------------------------

figure; hold on;
title('Training Loss')
xlabel('#Epoch')
ylabel('Training Loss')
plot(1:length(loss.train),loss.train,DisplayName='Training Loss');
legend

% Visualize the Output Space. ---------------------------------------------

% Compute output samples.
ys = nn.evaluate(xs);
% Obtain predicted classes.
[~,k] = max(ys,[],1);

% Plot everything.
figure; hold on;
plotPoints(xs(:,k==1),1:2,DisplayName='Predictions (Class 1)')
plotPoints(xs(:,k==2),1:2,DisplayName='Predictions (Class 2)')

% Set the output.
res = true;

end


% Auxiliary functions -----------------------------------------------------

function [xs,ts] = aux_generateTwoMoons(N,noise)
    % Generate a two moons dataset.
    % N: number of samples per moon
    % noise: standard deviation of Gaussian noise
    
    % Generate upper moon.
    t1 = linspace(0,pi,N);
    x1 = [cos(t1); sin(t1)] + noise*randn([2 N]);
    
    % Generate lower moon (shifted and flipped).
    t2 = linspace(0,pi,N);
    x2 = [1 - cos(t2); 0.5 - sin(t2)] + noise*randn([2 N]);
    
    % Combine the samples.
    xs = [x1 x2];
    ts = repelem(eye(2),1,N);
end

function [nn,layers,options] = aux_generateNeuralNetwork(options,n0,nk,nK,K,actfun)
    % Generate a random neural network.
    % - n0: number of input dimensions.
    % - nk: number of hidden dimensions.
    % - nK: number of output dimensions.
    % - K: number of layers.
    % - actfun: type of activation function, 
    %       i.e., {'relu','tanh',...}.

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
        layers{end+1} = nnLinearLayer(zeros(nout,nin),zeros(nout,1));
        if i < length(nks)-2
            % Create an activation layer.
            actl = nnActivationLayer.instantiateFromString(actfun);
            % Append the activation layer.
            layers{end+1} = actl;
        end
    end

    % Initialize the neural network.
    nn = neuralNetwork(layers);
    nn.setInputSize([n0 1]);
    % Initialize the weights and bias.
    nn.initWeights('glorot');
end

% ------------------------------ END OF CODE ------------------------------
