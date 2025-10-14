function res = test_nn_neuralNetwork_normWeights()
% test_nn_neuralNetwork_normWeights - tests the normWeights function 
%    
%
% Syntax:
%    res = test_nn_neuralNetwork_normWeights()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
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

% We construct a monotonic neural network and check if it is monotonic.

n0 = 1; % number of input dimensions.
nK = 1; % number of output dimensions.
nk = randi(100,1); % number of hidden dimensions.
K = randi(10,1); % number of layers.
actfun = 'groupSort'; % type of activation function, i.e., {'relu','tanh','groupSort'}.
lambda = randi(100,1); % Lipschitz constant
[nn,~] = aux_generateMonotonicNeuralNetwork(n0,nk,nK,K,actfun,lambda);

% Normalize the weights.
nn.normWeights();

% Sample random inputs.
N = 1000;
xs = linspace(-1,1,N);
% Compute outputs.
options.nn.train.backprop = true;
ys = nn.evaluate(xs,options);
[~,idx] = sort(ys,'ascend');

assert(all(idx == 1:N,'all'));

res = 1;

figure; hold on;
title('Monotonic Neural Network')
xlabel('Input')
ylabel('Output')
plot(xs,ys,'DisplayName','Neural Network');
legend

end


% Auxiliary functions -----------------------------------------------------

function [nn,layers] = aux_generateMonotonicNeuralNetwork(n0,nk,nK,K,actfun,lambda)
    % Generate a random neural network.
    % - n0: number of input dimensions.
    % - nk: number of hidden dimensions.
    % - nK: number of output dimensions.
    % - K: number of layers.
    % - actfun: type of activation function, 
    %       i.e., {'relu','tanh','groupSort'}.
    % - lambda: Lipschitz constant
    
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
        % Create a linear layer with random weights.
        layers{end+1} = nnLipConstrLinearLayer(ones(nout,nin),ones(nout,1),lambda^-(1/K));
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

end

% ------------------------------ END OF CODE ------------------------------
