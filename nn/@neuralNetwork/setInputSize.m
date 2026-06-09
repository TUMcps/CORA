function outputSize = setInputSize(obj, varargin)
% setInputSize - propagate inputSize through the network and store the
% inputSize for each layer. This is necessary to propagate images through a
% network.
%
% Syntax:
%    outputSize = setInputSize(obj, inputSize)
%
% Inputs:
%    inputSize - column vector, with sizes of each dimension
%    verbose: bool if information should be displayed
%    idxLayer: indices to layer, for which we set the input size
%
% Outputs:
%    outputSize - output size of the neural network
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: NeuralNetwork

% Authors:       Lukas Koller, Tobias Ladner
% Written:       10-June-2022
% Last update:   17-January-2023 (TL, polish)
% Last revision: 17-August-2022

% ------------------------------ BEGIN CODE -------------------------------

% Set default values.
[inputSize,verbose,idxLayer] = setDefaultValues({[],false,1:length(obj.layers)},varargin);

if isempty(inputSize)
    if ~isempty(obj.neurons_in)
        % Set default input size.
        inputSize = [obj.neurons_in 1];
    else
        throw(CORAerror('CORA:specialError', ...
            'Please provide an input size. Unable to determine it from network weights.'));
    end
end

% compute in-/out sizes of all layers
if verbose
    disp("Computing in-/out sizes of all layers...")
end
obj.neurons_in = prod(inputSize);
for i = idxLayer
    % iterate through all layers
    layer_i = obj.layers{i};
    outputSize = layer_i.computeSizes(inputSize);
    if verbose
        fprintf(" (%d)\t %s\n", i, layer_i.getLayerInfo())
    end
    inputSize = outputSize;
end
obj.neurons_out = prod(inputSize);
end

% ------------------------------ END OF CODE ------------------------------
