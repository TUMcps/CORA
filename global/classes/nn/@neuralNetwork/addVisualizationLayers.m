function nn = addVisualizationLayers(obj, varargin)
% addVisualizationLayers - insert visualization layers
%
% Syntax:
%    nn = addVisualizationLayers(obj, verbose)
%
% Inputs:
%    obj - neuralNetwork
%    neuronIds - 2 column array, ids of neurons to visualize
%    verbose - bool if information should be displayed
%
% Outputs:
%    nn - new neural network with visualization layers inserted
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: NeuralNetwork

% Authors:       Lukas Koller, Tobias Ladner
% Written:       23-June-2022
% Last update:   17-July-2023 (TL, parse input, clean up)
% Last revision: 17-August-2022

% ------------------------------ BEGIN CODE -------------------------------

% parse input
n = length(obj.layers);
[neuronIds, verbose] = setDefaultValues({[1 2],false},varargin);

inputArgsCheck({ ...
    {obj, 'att', 'neuralNetwork'}; ...
    {neuronIds, 'att', 'numeric', {'integer', 'positive'}}; ...
    {verbose, 'att', 'logical'}
})

% check dimensions
if size(neuronIds, 2) ~= 2
    throw(CORAerror("CORA:specialError", 'Specified neuronIds need to have two columns.'))
end
if size(neuronIds, 1) == 1
    % expand to number of layers
    neuronIds = repmat(neuronIds,n+1,1);
end
if size(neuronIds, 1) ~= n+1
    throw(CORAerror("CORA:specialError", 'Specified neuronIds need to match number of layers + 1.'))
end

% inject visualization layer
lastLayerName = "Input";
newLayers = {};
for i = 0:n
    if i > 0
        newLayers{2*i} = obj.layers{i};
    end

    if i < n
        nextLayerName = obj.layers{i+1}.name;
    else
        nextLayerName = "Output";
    end

    visName = sprintf("Visualize: %s -> %s", lastLayerName, nextLayerName);
    newLayers{2*i+1} = nnVisualizationLayer(i+1, neuronIds(i+1, :), visName);
    lastLayerName = nextLayerName;
end
newLayers = reshape(newLayers, [], 1); % column vector

% create network
nn = neuralNetwork(newLayers);

% verbose
if verbose
    display(nn);
end

end

% ------------------------------ END OF CODE ------------------------------
