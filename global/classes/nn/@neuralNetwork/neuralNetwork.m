classdef neuralNetwork
% neuralNetwork - class that stores layer-based neural networks
%
% Syntax:
%    obj = neuralNetwork(layers)
%
% Inputs:
%    layers - cell-array storing layers of type nnLayer
%
% Outputs:
%    obj - generated object
%
% Example:
%    layers = cell(2, 1);
%    W = rand(1,10); b = rand(1,1);
%    layers{1} = nnLinearLayer(W, b);
%    layers{2} = nnSigmoidLayer();
%    nn = neuralNetwork(layers);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: NNLayer

% Author:       Tobias Ladner
% Written:      28-March-2022
% Last update:  30-March-2022
%               23-November-2022 (polish)
% Last revision:10-August-2022 (renamed)

%------------- BEGIN CODE --------------

properties
    layers = []

    neurons_in
    neurons_out
end
methods
    % constructor
    function obj = neuralNetwork(layers)
        % check inputs
        if nargin ~= 1
            throw(CORAerror('CORA:notEnoughInputArgs', 1));
        end

        if ~(isa(layers, 'cell') ...
                && all(arrayfun(@(layer) isa(layer{1}, 'nnLayer'), layers)))
            throw(CORAerror('CORA:wrongInputInConstructor', ...
               'First argument should be a cell array of type nnLayer.'));
        end
        obj.layers = reshape(layers, [], 1);

        for i = 1:size(obj.layers, 1)
            [nin, ~] = obj.layers{i}.getNumNeurons();
            if ~isempty(nin)
                obj.neurons_in = nin;
                break;
            end
        end

        for i = size(obj.layers, 1):-1:1
            [~, nout] = obj.layers{i}.getNumNeurons();
            if ~isempty(nout)
                obj.neurons_out = nout;
                break;
            end
        end
    end

    % methods
    r = evaluate(obj, input, evParams)
end

methods (Static)
    obj = generateRandom(varargin)

    % conversions methods
    obj = readNetwork(file_path)
    obj = readONNXNetwork(file_path, verbose, inputDataFormats, outputDataFormats)
    obj = readNNetNetwork(file_path)
    obj = readYMLNetwork(file_path)
    obj = readSherlockNetwork(file_path)
    obj = convertDLToolboxNetwork(dltoolbox_layers, verbose)

    obj = getFromOldNeuralNetwork(nn_old)
end
end

%------------- END OF CODE --------------