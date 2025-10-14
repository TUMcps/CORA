function [res] = test_nn_neuralNetwork_getInputNeuronOrder()
% test_nn_neuralNetwork_getInputNeuronOrder - tests the 
%    getInputNeuronOrder function 
%    
%
% Syntax:
%    res = test_nn_neuralNetwork_getInputNeuronOrder()
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

% Authors:       Tobias Ladner
% Written:       08-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


seed = 1;
rng(seed);

% create simple network
W1 = rand(10, 16) *2 -1;
b1 = rand(10, 1);

W2 = rand(3, 10) *2 -1;
b2 = rand(3, 1); 

nn = neuralNetwork({ ...
    nnLinearLayer(W1, b1);
    nnSigmoidLayer();
    nnLinearLayer(W2, b2);
    nnSigmoidLayer();
});

x = rand(16,1);

% test each method ---

% in-order
neuronOrder = nn.getInputNeuronOrder('in-order',x);
assert(isequal(neuronOrder,1:16))

% sensitivity
neuronOrder = nn.getInputNeuronOrder('sensitivity',x);
[~,idx] = sort(mean(abs(nn.calcSensitivity(x)),1));
assert(isequal(neuronOrder,idx));

% snake
neuronOrder = nn.getInputNeuronOrder('snake',x,[4,4,1]);
assert(isequal(neuronOrder,[1, 5, 9, 13, 14, 15, 16, 12, 8, 4, 3, 2, 6, 10, 11, 7]))


% test completed
res = true;


end

% ------------------------------ END OF CODE ------------------------------
