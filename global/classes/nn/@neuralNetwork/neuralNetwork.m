classdef neuralNetwork < handle
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
% See also: nnLayer

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   30-March-2022
%                23-November-2022 (polish)
%                17-January-2023 (try to set input size)
% Last revision: 10-August-2022 (renamed)

% ------------------------------ BEGIN CODE -------------------------------

properties
    layers = []

    neurons_in
    neurons_out
end
methods
    % constructor
    function obj = neuralNetwork(layers)
        % check inputs
        if nargin < 1
            layers = {};
        elseif nargin > 1
            throw(CORAerror('CORA:tooManyInputArgs', 1));
        end

        if ~(isa(layers, 'cell') ...
                && all(arrayfun(@(layer) isa(layer{1}, 'nnLayer'), layers)))
            throw(CORAerror('CORA:wrongInputInConstructor', ...
               'First argument should be a cell array of type nnLayer.'));
        end
        obj.layers = reshape(layers, [], 1);

        % simple neurons_in and _out computation
        for i = 1:length(obj.layers)
            [nin, ~] = obj.layers{i}.getNumNeurons();
            if ~isempty(nin)
                obj.neurons_in = nin;
                break;
            end
        end
        for i = length(obj.layers):-1:1
            [~, nout] = obj.layers{i}.getNumNeurons();
            if ~isempty(nout)
                obj.neurons_out = nout;
                break;
            end
        end

        try            
            % to automatically determine input sizes of all layers
            obj.setInputSize();
        end
    end

    % evaluate ------------------------------------------------------------

    r = evaluate(obj, input, varargin)
    S = calcSensitivity(obj, x, varargin)
    
    % refine --------------------------------------------------------------
    
    refine(obj, max_order, type, method, x, verbose, force_bounds, gamma)
    refinable_layers = getRefinableLayers(obj)

    % reset ---------------------------------------------------------------

    reset(obj)
    resetApproxOrder(obj)
    resetBounds(obj)

    % verify --------------------------------------------------------------
    
    [res, x] = verify(nn, X0, spec, varargin)

    % Auxiliary functions -------------------------------------------------

    pattern = getOrderPattern(obj)
    numNeurons = getNumNeurons(obj)
    nn_normal = getNormalForm(obj)

    function l = length(obj)
        % returns the number of layers
        l = length(obj.layers);
    end
end

% Static functions --------------------------------------------------------

methods (Static)
    obj = generateRandom(varargin)

    % read & convert ------------------------------------------------------
    
    obj = readNetwork(file_path)
    obj = readONNXNetwork(file_path, varargin)
    obj = readNNetNetwork(file_path)
    obj = readYMLNetwork(file_path)
    obj = readSherlockNetwork(file_path)
    obj = readGNNnetwork(file_path, varargin)
    
    % convert
    obj = convertDLToolboxNetwork(dltoolbox_layers, verbose)
    obj = getFromCellArray(W, b, actFun)
end

% internal functions ------------------------------------------------------

methods (Access=protected)
    propagateBounds(obj, i, evParams)
end

end

% ------------------------------ END OF CODE ------------------------------
