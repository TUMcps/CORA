function res = test_nn_neuralNetwork_enumerateLayers()
% test_nn_neuralNetwork_enumerateLayers - unit test function for 
%     neuralNetwork/enumerateLayers
%
% Syntax:
%    res = test_nn_neuralNetwork_enumerateLayers()
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
% See also: neuralNetwork

% Authors:       Lukas Koller
% Written:       05-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Set result.
res = true;

% Test case 1: Sequence of layers. ----------------------------------------
% Specify the number of layers.
numLinLayers = 10;
% Generate a random neural network.
nn = neuralNetwork.generateRandom(NrLayers=numLinLayers);
% Obtain the number of layers.
K = length(nn.layers);
% Enumerate the layers.
[layersEnum,ancIdx,predIdx,succIdx] = nn.enumerateLayers();
% The neural network is a squence of layers; therefore, the enumerated
% layers are equal to the stored layers cell array.
assert(isequal(layersEnum,nn.layers'));
assert(all(ancIdx == 1:K));
assert(all(predIdx == 0:(K-1)));
assert(all(succIdx == 2:(K+1)));

% Test case 2: Composite layers. ------------------------------------------
nn = neuralNetwork({ ...
    nnLinearLayer(0,0); ...
    nnReLULayer; ...
    nnCompositeLayer({ ...
        {nnLinearLayer(0,0) nnReLULayer nnLinearLayer(0,0) nnReLULayer}; ...
        {nnLinearLayer(0,0) nnReLULayer}},'add');
    nnLinearLayer(0,0); ...
    nnReLULayer; ...
    nnLinearLayer(0,0)
});
% Enumerate the layers.
[layersEnum,ancIdx,predIdx,succIdx] = nn.enumerateLayers();

assert(isequal(layersEnum, ...
    [nn.layers(1:2)' ...
    nn.layers{3}.layers{1} ...
    nn.layers{3}.layers{2} ...
    nn.layers(4:6)']));
assert(all(ancIdx == [1:2 repelem(3,6) 4:6]));
assert(all(predIdx == [0:5 2 7:10]));
assert(all(succIdx == [2:6 9 8 9:12]));

% ------------------------------ END OF CODE ------------------------------
