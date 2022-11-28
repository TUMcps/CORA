function obj = readONNXNetwork(file_path, verbose, inputDataFormats, outputDataFormats)
% readONNXNetwork - reads and converts a network saved in onnx format
%
% Syntax:
%    res = neuralNetwork.readONNXNetwork(file_path)
%
% Inputs:
%    file_path - path to file
%    verbose - bool if information should be displayed
%    inputDataFormats - dimensions of input e.g. 'BC' or 'BSSC'
%    outputDataFormats - see inputDataFormats
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralnetwork2cora

% Author:       Tobias Ladner
% Written:      30-March-2022
% Last update:  07-June-2022 (specify in- & outputDataFormats)
% Last revision:---

%------------- BEGIN CODE --------------

% validate parameters
if nargin < 2
    verbose = false;
end
if nargin < 3
    inputDataFormats = 'BC';
end
if nargin < 4
    outputDataFormats = 'BC';
end

try
    if verbose
        disp("Reading Network: Try #1 ...")
    end

    dltoolbox_net = importONNXNetwork(file_path, ...
        'InputDataFormats', inputDataFormats, 'OutputDataFormats', outputDataFormats);
    obj = neuralNetwork.convertDLToolboxNetwork(dltoolbox_net.Layers, verbose);
catch ex
    try
        if verbose
            disp("Reading Network: Try #2 ...")
        end

        dltoolbox_net = importONNXNetwork(input, 'OutputLayerType', 'regression');
        obj = neuralNetwork.convertDLToolboxNetwork(dltoolbox_net.Layers, verbose);
    catch
        try
            if verbose
                disp("Reading Network: Try #3 ...")
            end

            dltoolbox_net = importONNXNetwork(input, 'OutputLayerType', ...
                'regression', 'TargetNetwork', 'dlnetwork');
            % TODO: remove dependency on neuralNetworkOld
            nn_old = aux_constructFromParams(dltoolbox_net.Learnables.Value);
            obj = neuralNetwork.getFromOldNeuralNetwork(nn_old);
        catch
            try
                if verbose
                    disp("Reading Network: Try #4 ...")
                end

                L = importONNXLayers(input, 'OutputLayerType', ...
                    'regression', 'ImportWeights', true);
                obj = neuralNetwork.convertDLToolboxNetwork(L, verbose);
            catch
                rethrow(ex);
            end
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function NN = aux_constructFromParams(params)
% construct a neural network from the learnable parameters

% TODO: remove dependency on neuralNetworkOld

% extract weights and biases from the learnable parameters
W = {};
b = {};

for i = 1:length(params)
    if size(params{i}, 2) == 1
        b{end+1, 1} = double(extractdata(params{i}));
    else
        W{end+1, 1} = double(extractdata(params{i}))';
    end
end

% try to construct a neuralNetworkOld object
actFun = [repmat({'ReLU'}, [length(W) - 1, 1]); {'identity'}];

NN = neuralNetworkOld(W, b, actFun);

warning(['Could not determine activation functions,', ...
    ' so ReLUs are used as a default.']);
end

%------------- END OF CODE --------------