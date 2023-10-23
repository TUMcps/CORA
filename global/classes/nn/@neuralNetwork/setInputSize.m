function outputSize = setInputSize(obj, inputSize, verbose)
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

if nargin < 2
    if isempty(obj.neurons_in)
        error("Please provide an input size. Unable to determine it from network weights.")
    end
    inputSize = [obj.neurons_in, 1];
end
if nargin < 3
    verbose = false;
end

if verbose
    disp("Computing in-/out sizes of all layers...")
end

obj.neurons_in = prod(inputSize);
for i = 1:length(obj.layers)
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
