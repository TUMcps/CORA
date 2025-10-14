function res = test_nn_neuralNetwork_castWeights()
% test_nn_neuralNetwork_castWeights - tests the castWeights function 
%    
%
% Syntax:
%    res = test_nn_neuralNetwork_castWeights()
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

% Specify number of input and output dimensions.
n0 = 5;
nK = 7;
% Generate a random neural network.
nn = neuralNetwork.generateRandom( ...
    'NrInputs',n0, ...
    'NrOutputs',nK, ...
    'NrLayers',3, ...
    'NrHiddenNeurons',17 ...
);
% Append a softmax layer.
nn.layers{end+1} = nnSoftmaxLayer;

% Construct a sample input.
x = single(1);
nn.castWeights(x);

% Check that all weights are singles.
[layers,~] = nn.enumerateLayers();
for i=1:length(layers)
    names = layers{i}.getLearnableParamNames();
    for j=1:length(names)
        assert(isa(layers{i}.(names{j}),class(x)));
    end
end

% Construct a sample input.
x = double(1);
nn.castWeights(x);

% Check that all weights are doubles.
[layers,~] = nn.enumerateLayers();
for i=1:length(layers)
    names = layers{i}.getLearnableParamNames();
    for j=1:length(names)
        assert(isa(layers{i}.(names{j}),class(x)));
    end
end

res = 1;

end

% ------------------------------ END OF CODE ------------------------------
