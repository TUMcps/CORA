function [res] = test_nn_neuralNetwork_neuralNetwork()
% test_nn_neuralNetwork_neuralNetwork - tests constructor of neuralNetwork
%
% Syntax:
%    res = test_nn_neuralNetwork_neuralNetwork()
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
% Written:       23-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% check simple example
layers = cell(2, 1);
W = rand(4,3); b = rand(4,1);
layers{1} = nnLinearLayer(W, b);
layers{2} = nnTanhLayer();
nn = neuralNetwork(layers);
if ~(nn.neurons_in == 3 && nn.neurons_out == 4)
    res = false;
end

% check larger example
layers = cell(4, 1);
W1 = rand(10,2); b1 = rand(10, 1);
W2 = rand(3, 10); b2 = rand(3,1);
layers{1} = nnLinearLayer(W1, b1);
layers{2} = nnReLULayer();
layers{3} = nnLinearLayer(W2, b2);
layers{4} = nnSigmoidLayer();
nn = neuralNetwork(layers);
if ~(nn.neurons_in == 2 && nn.neurons_out == 3)
    res = false;
end

try
    nn = neuralNetwork(nnSigmoidLayer()); % should throw an error.
    res = false;
catch
end


end

% ------------------------------ END OF CODE ------------------------------
