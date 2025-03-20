function res = test_nn_random()
% test_nn_random - tests a random network using a numeric input
%
% Syntax:
%    res = test_nn_random()
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
% See also: -

% Authors:       Tobias Ladner
% Written:       24-June-2022
% Last update:   28-November-2022 (name-value pair syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% check examples in syntax
nn = neuralNetwork.generateRandom();

% nrOfInputs
nrOfInputs = 2;
nn = neuralNetwork.generateRandom('NrInputs',nrOfInputs);
assert(nn.neurons_in == nrOfInputs);

% nrOfOutputs
nrOfOutputs = 3;
nn = neuralNetwork.generateRandom('NrInputs',nrOfInputs,...
  'NrOutputs',nrOfOutputs);
assert(nn.neurons_in == nrOfInputs & nn.neurons_out == nrOfOutputs);

% activation function
actFun = 'ReLU';
nn = neuralNetwork.generateRandom('NrInputs',nrOfInputs,...
  'NrOutputs',nrOfOutputs,'ActivationFun',actFun);
assert(nn.neurons_in == nrOfInputs & nn.neurons_out == nrOfOutputs ...
    & isa(nn.layers{2}, 'nnReLULayer'));

% nrOfLayers
nrLayers = 3;
nn = neuralNetwork.generateRandom('NrInputs',nrOfInputs,...
  'NrOutputs',nrOfOutputs,'ActivationFun',actFun,'NrLayers',nrLayers);
assert(nn.neurons_in == nrOfInputs & nn.neurons_out == nrOfOutputs ...
    & isa(nn.layers{2}, 'nnReLULayer') * length(nn.layers) == 2*nrLayers);

% check other example
actFun = 'tanh';
nrOfInputs = 4;
nn = neuralNetwork.generateRandom('NrInputs',nrOfInputs, ...
    'ActivationFun',actFun);
assert(nn.neurons_in == nrOfInputs & isa(nn.layers{2}, 'nnTanhLayer'));

end

% ------------------------------ END OF CODE ------------------------------
