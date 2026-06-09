function res = test_nn_neuralNetwork_generateFromString()
% test_nn_neuralNetwork_generateFromString - tests the construction of a
%    neural network from a string specifying the layers
%
% Syntax:
%    res = test_nn_neuralNetwork_generateFromString()
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
% See also: neuralNetwork/generateFromString

% Authors:       Lukas Koller
% Written:       17-April-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Specify the layers of a small conv. network for MNIST.
layersstring = [
    "CONV 5 4*4+2"
    "CONV 10 4*4+1"
    "FC 100"
];
inSz = [28 28 1];
nK = 10;

% Generate the neural network.
nn = neuralNetwork.generateFromString(layersstring,inSz,nK);

% Check input/output dimensions.
assert(nn.neurons_in == prod(inSz));
assert(nn.neurons_out == nK);

% Check the composition of the layers: for each of the 3 specified layers
% we expect a linear/convolutional layer followed by an activation layer;
% additionally, there is a final linear output layer.
assert(length(nn.layers) == 2*3 + 1);
assert(isa(nn.layers{1},'nnConv2DLayer'));
assert(isa(nn.layers{2},'nnReLULayer'));
assert(isa(nn.layers{3},'nnConv2DLayer'));
assert(isa(nn.layers{4},'nnReLULayer'));
assert(isa(nn.layers{5},'nnLinearLayer'));
assert(isa(nn.layers{6},'nnReLULayer'));
assert(isa(nn.layers{7},'nnLinearLayer'));

% Generate a feed-forward network with tanh activations and batch norm.
layersstring = [
    "FC 20"
    "FC 15"
];
nn = neuralNetwork.generateFromString(layersstring,[4 1],3, ...
    'tanh',true);

% Check input/output dimensions.
assert(nn.neurons_in == 4);
assert(nn.neurons_out == 3);

% Expect a batch norm + activation layer after each FC layer.
assert(length(nn.layers) == 3*2 + 1);
assert(isa(nn.layers{1},'nnLinearLayer'));
assert(isa(nn.layers{2},'nnBatchNormLayer'));
assert(isa(nn.layers{3},'nnTanhLayer'));
assert(isa(nn.layers{4},'nnLinearLayer'));
assert(isa(nn.layers{5},'nnBatchNormLayer'));
assert(isa(nn.layers{6},'nnTanhLayer'));
assert(isa(nn.layers{7},'nnLinearLayer'));

% Test that a string containing newlines (and blank lines) is parsed
% correctly and empty lines are ignored.
layersstring = sprintf('FC 8\n\nFC 4\n');
nn = neuralNetwork.generateFromString(layersstring,[2 1],2);
assert(nn.neurons_in == 2);
assert(nn.neurons_out == 2);
assert(length(nn.layers) == 2*2 + 1);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
