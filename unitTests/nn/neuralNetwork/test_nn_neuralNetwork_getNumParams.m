function res = test_nn_neuralNetwork_getNumParams()
% test_nn_neuralNetwork_getNumParams - unit test function for 
%     neuralNetwork/getNumParams
%
% Syntax:
%    res = test_nn_neuralNetwork_getNumParams()
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
% See also: neuralNetwork/evaluate

% Authors:       Lukas Koller
% Written:       15-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Test empty case.
nn = neuralNetwork();
% Compute the number of parameters.
numParams = nn.getNumParams();
assert(isscalar(numParams) && numParams == 0);

% Activation layers do not have any learnable parameters.
nn = neuralNetwork({nnSigmoidLayer()});
% Compute the number of parameters.
numParams = nn.getNumParams();
assert(isscalar(numParams) && numParams == 0);

% Activation layers do not have any learnable parameters.
nn = neuralNetwork({nnReLULayer()});
% Compute the number of parameters.
numParams = nn.getNumParams();
assert(isscalar(numParams) && numParams == 0);

% Activation layers do not have any learnable parameters.
nn = neuralNetwork({nnTanhLayer()});
% Compute the number of parameters.
numParams = nn.getNumParams();
assert(isscalar(numParams) && numParams == 0);

% Test single layer.
n0 = randi(100); % Random number of input neurons.
nK = randi(100); % Random number of output neurons.
% Initialize a single linear layer neural network.
nn = neuralNetwork({nnLinearLayer(zeros([nK n0]),zeros([nK 1]))});
% Compute the number of parameters.
numParams = nn.getNumParams();
assert(numParams == (n0*nK) + nK);

% Initialize a single linear layer without learnable weight.
nn = neuralNetwork({nnLinearLayer(zeros([nK n0]),zeros([nK 1]),'name',false)});
% Compute the number of parameters.
numParams = nn.getNumParams();
assert(numParams == 0);

% Generate a random larger neural network.
K = randi(20); % Random number of layers.
nk = randi(100); % Random number of hidden neurons.
nn = neuralNetwork.generateRandom('NrInputs',n0,'NrOutputs',nK, ...
    'ActivationFun','relu','NrLayers',K,'NrHiddenNeurons',nk);
% Compute the number of parameters.
numParams = nn.getNumParams();
if K == 1
    assert(numParams == (n0+1)*nK);
else
    assert(numParams == (n0+1)*nk + (K-2)*(nk+1)*nk + (nk+1)*nK);
end

% Initialize a single convolutional layer neural network.
kh = randi(10); % Random kernel height.
kw = randi(10); % Random kernel width.
cin = randi(10); % Random number of input channels.
cout = randi(10); % Random number of output channels.
nn = neuralNetwork({nnConv2DLayer(zeros([kh kw cin cout]),zeros([cout 1]))});
% Compute the number of parameters.
numParams = nn.getNumParams();
assert(numParams == (kh*kw*cin+1)*cout);

% gather results
res = true;


% ------------------------------ END OF CODE ------------------------------
