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
%    W = [ -0.928, 1.244, -0.589, 1.491, -0.806, -1.113, 0.957, -1.363, 0.681, 0.438 ]; b = -0.143;
%    layers{1} = nnLinearLayer(W, b);
%    layers{2} = nnSigmoidLayer();
%    nn = neuralNetwork(layers);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnLayer

% Authors:       Tobias Ladner, Lukas Koller
% Written:       28-March-2022
% Last update:   30-March-2022
%                23-November-2022 (polish)
% Last revision: 10-August-2022 (renamed)
%                02-August-2023 (LK, headers)

% ------------------------------ BEGIN CODE -------------------------------

properties
    layers = []

    neurons_in
    neurons_out

    reductionRate = 1;
end
methods
    % constructor
    function obj = neuralNetwork(layers)
        % check inputs
        narginchk(0,1);
        if nargin < 1
            layers = {};
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

    r = evaluate_(obj, input, options, idxLayer)
    [S,y] = calcSensitivity(obj, x, varargin)
    
    % refine --------------------------------------------------------------
    
    refine(obj, max_order, type, method, x, verbose, force_bounds, gamma)
    refinable_layers = getRefinableLayers(obj)

    % reset ---------------------------------------------------------------

    reset(obj)
    resetApproxOrder(obj)
    resetBounds(obj)
    resetGNN(obj)

    % verify --------------------------------------------------------------
    
    [res, x_, y_] = verify(nn, x, r, A, b, safeSet, varargin)

    % reduce --------------------------------------------------------------
    
    [nn_red, S] = computeReducedNetwork(obj, S, varargin)

    % training ------------------------------------------------------------

    initWeights(nn, varargin)

    castWeights(nn, varargin)

    normWeights(nn, options, varargin)

    [loss,trainTime] = train(nn, trainX, trainY, valX, valY, options, verbose)

    [grad_in] = backprop(nn, grad_out, options, varargin)

    [l, u] = backpropIntervalBatch(nn, l, u, options, varargin)

    [numGen,nn] = prepareForZonoBatchEval(nn, x, varargin)
    [c, G] = evaluateZonotopeBatch(nn, c, G, varargin)
    [c, G] = evaluateZonotopeBatch_(nn, c, G, options, idxLayer)
    [gc, gG] = backpropZonotopeBatch(nn, gc, gG, varargin)
    [gc, gG] = backpropZonotopeBatch_(nn, gc, gG, options, idxLayer, updateWeights)

    % explain -------------------------------------------------------------

    [idxFreedFeats,featOrder,timesPerFeat] = explain(nn, x, target, epsilon, varargin)

    % convert & export ----------------------------------------------------

    nn_dlt = convertToDLToolboxNetwork(nn)
    res = exportONNXNetwork(nn,file_path,varargin)
    res = exportNetworkAsCellArray(nn,file_path)
    jsonstr = exportAsJSON(obj,varargin);
    nnStruct = exportAsStruct(obj);

    % Auxiliary functions -------------------------------------------------

    pattern = getOrderPattern(obj)
    numNeurons = getNumNeurons(obj)
    nn_normal = getNormalForm(obj)
    [layersEnum,ancIdx,predIdx,succIdx] = enumerateLayers(obj)
    neuronOrder = getInputNeuronOrder(obj,method,x,inputSize)
    gnn_red = reduceGNNForNode(obj,G,n0)

    function l = length(obj)
        % returns the number of layers
        l = length(obj.layers);
    end
end

% Static functions --------------------------------------------------------

methods (Static)
    obj = generateRandom(varargin)

    % read & convert ------------------------------------------------------
    
    % read
    obj = readNetwork(file_path)
    obj = readONNXNetwork(file_path, varargin)
    obj = readNNetNetwork(file_path)
    obj = readYMLNetwork(file_path)
    obj = readSherlockNetwork(file_path)
    obj = readGNNnetwork(file_path, varargin)
    data = readGNNdata(file_path, varargin)
    obj = readJSONNetwork(file_path)

    % import
    obj = importFromJSON(jsonstr);
    obj = importFromStruct(nnStruct);
    
    % convert
    obj = convertDLToolboxNetwork(dltoolbox_layers, verbose)
    obj = getFromCellArray(W, b, actFun)
end

% internal functions ------------------------------------------------------

methods (Access=protected)
    propagateBounds(obj, i, options)
end

end

% ------------------------------ END OF CODE ------------------------------
